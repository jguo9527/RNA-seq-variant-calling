import os
import subprocess

def mark_duplicates(input_bam, output_bam, metrics_file, picard_jar):
    """
    Run Picard MarkDuplicates on a BAM file.
    
    :param input_bam: Path to the input BAM file.
    :param output_bam: Path to the output BAM file (with duplicates marked).
    :param metrics_file: Path to the file where duplication metrics will be saved.
    :param picard_jar: Path to the Picard JAR file.
    """
    command = [
        "java", "-jar", picard_jar, "MarkDuplicates",
        "I=" + input_bam,  # Input BAM file
        "O=" + output_bam,  # Output BAM file
        "M=" + metrics_file,  # Metrics file for duplication statistics
        "REMOVE_DUPLICATES=true"  # Remove duplicates (set to false if you want to keep them)
    ]
    
    try:
        # Run the Picard command
        print(f"Running command: {' '.join(command)}")
        subprocess.run(command, check=True)
        print(f"Marked duplicates for {input_bam}, output saved to {output_bam}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {input_bam}: {e}")

def process_bam_files(input_folder, output_folder, picard_jar):
    """
    Process all BAM files in a folder and run MarkDuplicates on each.
    
    :param input_folder: Folder containing BAM files to process.
    :param output_folder: Folder to save the output BAM files with marked duplicates.
    :param picard_jar: Path to the Picard JAR file.
    """
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over each file in the folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".bam"):
            input_bam = os.path.join(input_folder, file_name)
            output_bam = os.path.join(output_folder, f"marked_{file_name}")
            metrics_file = os.path.join(output_folder, f"{file_name}_metrics.txt")

            # Mark duplicates for the current BAM file
            print(f"Processing: {file_name}")
            mark_duplicates(input_bam, output_bam, metrics_file, picard_jar)

if __name__ == "__main__":
    # Path to the folder containing the BAM files
    input_folder = "/Volumes/KIOXIA/Hemp/alignment/RG_added"  # Replace with the actual folder path
    output_folder = "/Volumes/KIOXIA/Hemp/alignment/duplicates_marked"  # Replace with the desired output folder path
    
    # Path to the Picard JAR file
    picard_jar = "/Volumes/KIOXIA/Hemp/alignment/picard.jar"  # Replace with the actual Picard JAR file path
    
    # Run the process
    process_bam_files(input_folder, output_folder, picard_jar)