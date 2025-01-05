#!/usr/bin/perl -w

###########################################################################################################
# Process demultiplexed reads, map, and extract allele counts at target positions
# 1) Align processed sequences using BWA and output SAM file
# 2) Edit SAM to solve soft clipped bases at the end of a fragment (call noEndClip2.py) (not required)
# 3) Convert SAM to BAM and then pileup
# 4) Extract depth from mpileup at target positions (call allelecount2.py)
###########################################################################################################
##########################################################################################
# Requirement 1: BWA aligner (version and parameters need to be decided)
# Requirement 2: SAMtools (Li et al., 2009)
##########################################################################################
use strict;
#use warnings;
no warnings 'uninitialized';
use Getopt::Long qw(GetOptions);
use Getopt::Long qw(:config no_ignore_case);
#use List::Util qw( min max );
use Parallel::ForkManager;

######################
# The help function
######################
my $help = $ARGV[0];
my ($bwa,$acpy,$samtools,$Reference,$outPos,$phred_Q,$map_q,$F,$f,$threads,$sam_add);

my $H = "\n###########################################################\n"
	."BWA alignment pipeline\n"
	."###########################################################\n"
	."Version: 1\n"
	."Align with BWA-mem and process with SAMtools\n"
	."Options:\n"
	."-bw: Path to BWA executable file. String."
	."-acpy: Path to python program, allelecount2."
	."-st: Path to SAMtools executable file. String."
	."-ref: Reference FASTA file, either Mock Reference or true reference."
	."-outPos: target SNP position file. *.txt\n"
	."-Q: Phred score base call quality. Numeric. "
	."-q: Alignment quality. Numeric. "
	."-F: SAMtools flags controlled by CAPS F (i.e. flag reads to exclude). Numeric."
	."-f: SAMtools flags controlled by small f (i.e. required flags to keep reads). Numeric. "
	."-t: Number of independent threads used. Numeric. "
	."-opt: If desired, any additional options for SAMtools view. String within “quotes”. Default: 0 (nothing)\n\n";

if (! defined $help or $help =~ "h" or $help =~ "H")  {
	print "$H";
	goto FINAL;
}

#################################
# Setting the parameter values
#################################
$bwa = '~/bwa-mem2.2.1/bwa-mem2';
$acpy = '~/allelecount/allelecount2.py';
$samtools = '~/samtools/samtools';
$Reference = 'Genome.fa';
$phred_Q = 0;
$map_q = 0;
$F = 4;
#$f = 4095;
$threads = 8;
$sam_add = 0;

GetOptions(
'bw=s' => \$bwa,	          	# string
'acpy=s' => \$acpy,
'st=s' => \$samtools,          	# string
'ref=s' => \$Reference,         # file
'outPos=s' => \$outPos,					# file
'Q=s' => \$phred_Q,             # numeric
'q=s' => \$map_q,               # numeric
'F=s' => \$F,                	# numeric
'f=s' => \$f,               	# numeric
't=s' => \$threads,             # numeric
'opt=s' => \$sam_add,           # string
) or die "$H\n";

print "\n#################################\n# Fastq to read depth, v1\n#################################\n";
my $pm = new Parallel::ForkManager($threads);
my $sttime = time;

print "\n Parameters used: \n
   bwa version = $bwa \n
   acpy version = $acpy \n
   samtools version = $samtools \n
   reference used = $Reference \n
   target positions = $outPos \n
   phred quality Q = $phred_Q \n
   mapping quality q = $map_q \n
   samtools F flag = $F \n
   samtools f flag = $f \n
   threads used = $threads \n
";

# =========================Prepare file names and space========================= #
# Create alignment directory
my $dir_align = "alignments";
unless(-e $dir_align, or mkdir $dir_align) {die "Directory $dir_align cannot be created.\n";}

# Create read depth (at outPos) directory
my $dir_allele = "alleledepth";
unless(-e $dir_allele, or mkdir $dir_allele) {die "Directory $dir_allele cannot be created.\n";}

# Get sample names
my @files = glob "*.bam";
my @samples =();
foreach my $sample (@files){
	push @samples, substr($sample, 0, index($sample, '__Aligned.sortedByCoord.out.bam'))
}
# write sample names to a text file
open my $sampleList, ">", "Sample_list.txt" or die "Can't initialize Sample_list.txt file\n";
foreach my $sample (@samples){
	print $sampleList "$sample\n";
}
# import sample names to files
chomp (@samples);
#print @samples;

# Indexing the reference FASTA sample
print "\nIndexing the reference genome FASTA file ...";
system ( "samtools faidx $Reference" );
print "\nDONE.\n\n";

# Pileup alleles and extract read dpeth at target postions
print "Extracting read depth at target positions ...\n";
foreach my $sample (@samples) {
	my $pid = $pm->start and next;
	my $bam_in = join ("$sample","__Aligned.sortedByCoord.out.bam");
	my $count_out = join (".","$sample","counts.txt");
	#system ( "samtools mpileup -Q $phred_Q -q $map_q -AB -x -f $Reference $input > $mpileup" );
	system ( "samtools mpileup -Q $phred_Q -q $map_q -AB -x -l $outPos -f $Reference $bam_in 2>/dev/null \\\
	| bcftools view -bvcg my-raw.bcf > my-var.bcf" );
	print "Extracting reads for sample $sample ...\n";
	$pm->finish;
}
$pm->wait_all_children;
print "DONE.\n\n";

# ================= get allele count for all samples ====================
print "Joining allele counts for all samples ...\n";
open my $POS, "<", "$outPos" or die "Can't load file $!";
open my $counts_out, ">Master_count.txt" or die "Can't load file $!";
# write scaffold and position to posHash
my %posHash;
while (<$POS>){
		my ($first, $second) =  split;
		$posHash{"$first $second"} = ();
}
close $POS;
# merge read depth data from all samples
foreach my $sample (@samples)  {
	open my $allele_count_file, "<", join(".","$sample","counts.txt") or die "Can't load file $!";

	my %genoHash;
	while ( <$allele_count_file> ){
		my ($first, $second, $count) =  split;
		$genoHash{"$first $second"} = $count;
	}

	foreach my $sc_pos ( keys %posHash ){
		push( @{ $posHash{$sc_pos} }, $genoHash{$sc_pos} ? $genoHash{$sc_pos} : "_,_,_,_" );
	}

	%genoHash = ();
	close $allele_count_file;
}
# print read depth data to master count file
foreach my $key (sort keys %posHash ) {
	print $counts_out join("\t", $key, @{ $posHash{$key} }),"\n" ;
}
close $counts_out;
print "DONE.\n";

system ( "mv *sorted.bam* ./alignments" );
system ( "mv *.counts.txt ./alleledepth" );
print "Elapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";

FINAL:
exit;
