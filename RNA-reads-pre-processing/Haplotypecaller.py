import os
import subprocess

def run_haplotype_caller(bam_file, reference_genome, output_dir, gatk_path, stand_call_conf=20):
    """
    Run HaplotypeCaller from GATK for a given BAM file with JVM options.

    Parameters:
    bam_file (str): Path to the BAM file.
    reference_genome (str): Path to the reference genome FASTA file.
    interval_list (str): Path to the interval list file.
    dbsnp_vcf (str): Path to the dbSNP VCF file.
    output_dir (str): Directory to store the output VCF files.
    gatk_path (str): Path to the GATK executable.
    stand_call_conf (int): The standard minimum confidence threshold for calling variants (default is 20).
    
    Returns:
    str: Output file path.
    """
    # Define the output VCF file path
    bam_filename = os.path.basename(bam_file)
    base_name = os.path.splitext(bam_filename)[0]
    output_vcf = os.path.join(output_dir, f"{base_name}.vcf.gz")

    # Construct the HaplotypeCaller command
    command = [
        "java",
        "-Xms6000m",  # Set initial heap size to 6GB
        "-XX:GCTimeLimit=50", 
        "-XX:GCHeapFreeLimit=10",
        "-jar", gatk_path,  # GATK executable path
        "HaplotypeCaller",  # Tool name
        "-R", reference_genome,  # Reference genome
        "-I", bam_file,  # Input BAM file
        "-O", output_vcf,  # Output VCF file
        "-dont-use-soft-clipped-bases",  # Option for base calling
        "--standard-min-confidence-threshold-for-calling", str(stand_call_conf),  # Calling confidence threshold
    ]

    # Run the command
    try:
        subprocess.run(command, check=True)
        print(f"HaplotypeCaller completed for {bam_file}. Output: {output_vcf}")
        return output_vcf
    except subprocess.CalledProcessError as e:
        print(f"Error running HaplotypeCaller for {bam_file}: {e}")
        return None

def process_bam_files(input_folder, reference_genome, output_dir, gatk_path, stand_call_conf=20):
    """
    Process all BAM files in a given folder, running HaplotypeCaller for each.

    Parameters:
    input_folder (str): Folder containing BAM files.
    reference_genome (str): Path to the reference genome FASTA file.
    interval_list (str): Path to the interval list file.
    dbsnp_vcf (str): Path to the dbSNP VCF file.
    output_dir (str): Directory to store output VCF files.
    gatk_path (str): Path to the GATK executable.
    stand_call_conf (int): The standard minimum confidence threshold for calling variants (default is 20).
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Loop through all BAM files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith(".bam"):
            bam_file = os.path.join(input_folder, filename)
            print(f"Processing {bam_file}...")
            run_haplotype_caller(bam_file, reference_genome, output_dir, gatk_path, stand_call_conf)

if __name__ == "__main__":
    # Define parameters
    input_folder = "/Volumes/KIOXIA/Hemp/alignment/Ncigar_split"  # Path to the folder containing BAM files
    reference_genome = "/Volumes/KIOXIA/Hemp/alignment/reference/GCF_029168945.1_ASM2916894v1_rna.fasta"  # Path to the reference genome FASTA file
    output_dir = "/Volumes/KIOXIA/Hemp/alignment/Vcf_called"  # Directory to store VCF output files
    gatk_path = "/Volumes/KIOXIA/Hemp/alignment/gatk-package-4.6.1.0-local.jar"  # Path to the GATK executable (e.g., /usr/local/bin/gatk)

    # Process BAM files in the input folder
    process_bam_files(input_folder, reference_genome, output_dir, gatk_path)
