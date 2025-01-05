import os
import subprocess

def split_ncigar_reads(input_bam, output_bam, gatk_jar, reference_fasta):
    """
    Run GATK SplitNCigarReads on a BAM file.
    
    :param input_bam: Path to the input BAM file.
    :param output_bam: Path to the output BAM file (with NCigar reads split).
    :param gatk_jar: Path to the GATK JAR file.
    :param reference_fasta: Path to the reference genome FASTA file.
    """
    command = [
        "java", "-jar", gatk_jar, "SplitNCigarReads",
        "-I", input_bam,  # Input BAM file
        "-O", output_bam,  # Output BAM file
        "-R", reference_fasta  # Reference genome in FASTA format
        # "--split-intervals",  # Optional, split intervals (you can use other options if needed)
        # "--skip-read-group-validation",  # Optional, skip read group validation (use as needed)
    ]
    
    try:
        # Run the GATK command
        print(f"Running command: {' '.join(command)}")
        subprocess.run(command, check=True)
        print(f"Successfully split NCigar reads for {input_bam}, output saved to {output_bam}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {input_bam}: {e}")

def process_bam_files(input_folder, output_folder, gatk_jar, reference_fasta):
    """
    Process all BAM files in a folder and run SplitNCigarReads on each.
    
    :param input_folder: Folder containing BAM files to process.
    :param output_folder: Folder to save the output BAM files with NCigar reads split.
    :param gatk_jar: Path to the GATK JAR file.
    :param reference_fasta: Path to the reference genome FASTA file.
    """
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over each file in the folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".bam"):
            input_bam = os.path.join(input_folder, file_name)
            output_bam = os.path.join(output_folder, f"split_{file_name}")

            # Run SplitNCigarReads for the current BAM file
            print(f"Processing: {file_name}")
            split_ncigar_reads(input_bam, output_bam, gatk_jar, reference_fasta)

if __name__ == "__main__":
    # Path to the folder containing the BAM files
    input_folder = "/Volumes/KIOXIA/Hemp/alignment/duplicates_marked"  # Replace with the actual folder path
    output_folder = "/Volumes/KIOXIA/Hemp/alignment/Ncigar_split"  # Replace with the desired output folder path
    
    # Path to the GATK JAR file
    gatk_jar = "/Volumes/KIOXIA/Hemp/alignment/gatk-package-4.6.1.0-local.jar"  # Replace with the actual path to the GATK JAR file
    
    # Path to the reference genome FASTA file
    reference_fasta = "/Volumes/KIOXIA/Hemp/alignment/reference/GCF_029168945.1_ASM2916894v1_genomic.fasta"  # Replace with the actual reference genome FASTA file path
    
    # Run the process
    process_bam_files(input_folder, output_folder, gatk_jar, reference_fasta)
