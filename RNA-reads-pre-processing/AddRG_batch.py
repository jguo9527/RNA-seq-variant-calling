import os
import subprocess

def add_read_group(input_bam, output_bam, read_group_template):
    """
    Adds or replaces the read group information in a BAM file using samtools addreplacerg.
    :param input_bam: Path to the input BAM file.
    :param output_bam: Path to the output BAM file.
    :param read_group_template: The template for the read group string, with a placeholder for sample name.
    """
    # Extract sample name from the BAM file name (remove the .bam extension)
    sample_name = os.path.splitext(os.path.basename(input_bam))[0]
    
    # Replace the placeholder with the actual sample name
    read_group = read_group_template.replace("{SM}", sample_name)

    # Construct the samtools command
    command = [
        "samtools", "addreplacerg",
        "-r", read_group,  # Add the read group information
        "-o", output_bam,  # Output file path
        input_bam  # Input BAM file path
    ]

    try:
        # Run the command
        print(f"Running command: {' '.join(command)}")
        subprocess.run(command, check=True)
        print(f"Read group added successfully: {output_bam}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {input_bam}: {e}")

def process_bam_files(folder_path, output_folder, read_group_template):
    """
    Process all BAM files in the folder and add read group information using samtools addreplacerg.
    :param folder_path: Path to the folder containing BAM files.
    :param output_folder: Path to the folder where output BAMs will be saved.
    :param read_group_template: The template for the read group information to be added.
    """
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over each file in the folder
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".bam"):
            input_bam = os.path.join(folder_path, file_name)
            output_bam = os.path.join(output_folder, f"rg_{file_name}")  # Prefix "rg_" for output files
            
            print(f"Processing: {file_name}")
            add_read_group(input_bam, output_bam, read_group_template)

if __name__ == "__main__":
    # Path to the folder containing the BAM files
    input_folder = "/Volumes/KIOXIA/Hemp/alignment/alignment-files"  # Replace with the actual folder path
    output_folder = "/Volumes/KIOXIA/Hemp/alignment/RG_added"  # Replace with the desired output folder path

    # The read group template with a placeholder for the sample name
    read_group_template = "@RG\tID:RG1\tSM:{SM}\tPL:Illumina\tLB:Library.fa"  # {SM} will be replaced by the BAM file name

    # Run the process
    process_bam_files(input_folder, output_folder, read_group_template)
