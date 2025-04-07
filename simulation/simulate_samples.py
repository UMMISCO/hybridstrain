""" 
Script to run the samples simulation
"""

import os
import subprocess
import shlex
import time
import logging

# configuring logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def nanosim_simulation(info_folder: str, out_dir: str, nanosim_model_and_prefix: str, replicates: int = 1, cpu: int = 20, parallel: int = 1, sample_range: tuple = None):
    """ 
    Generates simulated metagenomes using NanoSim, based on abundance tables in `info_folder`. The abundance tables should 
    each be named `<sample>.tsv` and a file named `genomes_list.tsv` should list the paths to genomes.

    Parameters:
    - info_folder (str): Path to the folder containing abundance tables and the genomes list (`genomes_list.tsv`).
    - out_dir (str): Path to the directory where NanoSim output files will be saved.
    - nanosim_model_and_prefix (str): Path to the NanoSim error model and its prefix (e.g., 'nanosim_models/metagenome_ERR3152364_Even/training').
    - replicates (int, optional): Number of replicate simulations to generate (default is 1).
    - cpu (int, optional): Number of CPU cores to use for parallel processing (default is 20).
    - parallel (int, optional): Number of parallel jobs to run with GNU parallel (default is 1).
    - sample_range (tuple, optional): Tuple specifying the range of samples to process in alphabetical order (default is None).

    Dependencies:
    - `NanoSim`: Ensure NanoSim is installed and `simulator.py` is accessible.
    - `pigz`: Fast parallel gzip compression is required for compressing FASTA output files.
    """

    samples_abundance = {}

    # collecting sample abundance files from info_folder, ignoring 'genomes_list.tsv'
    for filename in os.listdir(info_folder):
        if filename != "genomes_list.tsv" and filename.endswith(".tsv"):
            # extracting sample ID (filename without extension)
            sample_id = os.path.splitext(filename)[0]
            # adding full path of the abundance file to the dictionary
            samples_abundance[sample_id] = os.path.realpath(os.path.join(info_folder, filename))

    # sorting samples alphabetically and applying the sample range if provided
    sorted_samples = sorted(samples_abundance.keys())
    if sample_range:
        sorted_samples = sorted_samples[sample_range[0]:sample_range[1]]

    # logging the list of samples to process
    logging.info("Samples to process:")
    for sample_id in sorted_samples:
        logging.info(sample_id)
    
    time.sleep(3)

    # defining the path to the genome list file
    genomes_list = os.path.join(info_folder, "genomes_list.tsv")

    for i in range(replicates):
        commands = []
        # running NanoSim simulation for each sample in each replicate
        for sample in sorted_samples:
            abundance_path = samples_abundance[sample]
            logging.info(f"Processing sample: {sample}")
            
            # constructing output path prefix for each sample
            outdir_and_prefix = os.path.join(out_dir, str(i), sample)
            
            # defining the path for the combined reads file
            combined_reads = f"{outdir_and_prefix}_nanopore.fasta.gz"
            
            # check if the combined reads file already exists
            if os.path.exists(combined_reads):
                logging.info(f"Sample {sample} in replicate {i} already simulated. Skipping...")
                continue
            
            # assemblying NanoSim command with shell-safe arguments
            command = [
                "simulator.py", "metagenome", 
                "-gl", genomes_list, 
                "-a", abundance_path, 
                "-t", str(cpu),
                "-c", nanosim_model_and_prefix,
                "-o", outdir_and_prefix,
                "--seed", str(i)
            ]
            commands.append(command)
        
        # running all commands in parallel using GNU parallel
        if commands:
            logging.info(f"Running {len(commands)} NanoSim commands in parallel with {parallel} jobs.")
            parallel_command = f"parallel -j {parallel} ::: " + " ".join(shlex.quote(' '.join(cmd)) for cmd in commands)
            subprocess.run(parallel_command, shell=True, check=True)
        
        # post-processing for each sample
        for sample in sorted_samples:
            outdir_and_prefix = os.path.join(out_dir, str(i), sample)
            combined_reads_uncompressed = f"{outdir_and_prefix}_nanopore.fasta"
            aligned_reads = f"{outdir_and_prefix}_sample0_aligned_reads.fasta"
            unaligned_reads = f"{outdir_and_prefix}_sample0_unaligned_reads.fasta"
            error_profile = f"{outdir_and_prefix}_sample0_aligned_error_profile"

            # verifying and concatenating aligned and unaligned reads if both exist
            if os.path.exists(aligned_reads) and os.path.exists(unaligned_reads):
                logging.info(f"Both aligned and unaligned reads exist for sample {sample} in replicate {i}. Merging...")
                with open(combined_reads_uncompressed, 'w') as outfile:
                    for fname in [aligned_reads, unaligned_reads]:
                        with open(fname) as infile:
                            outfile.write(infile.read())
            
                # removing individual aligned and unaligned reads files after merging
                os.remove(aligned_reads)
                os.remove(unaligned_reads)
                # removing error profile file
                os.remove(error_profile)
                logging.info(f"Successfully merged and cleaned up reads for sample {sample} in replicate {i}")
            else:
                logging.warning(f"Warning: One or both FASTA files are missing for sample {sample} in replicate {i}")

            # compressing the combined FASTA file using pigz
            logging.info(f"Compressing combined reads file for sample {sample} in replicate {i}")
            compress_command = f"pigz {shlex.quote(combined_reads_uncompressed)}"
            subprocess.run(compress_command, shell=True, check=True)
            
def metagen_simulation(info_folder: str, out_dir: str, target_reads_number: int, ref_genomes: str,
                       read_length: int = 150, replicates: int = 1, 
                       cpu: int = 20):
    """ 
    Generates simulated metagenomes using the `metagen-simulator2` tool. This function reads abundance tables 
    stored in `info_folder`, where each table file is expected to be named `<sample>.tsv`.

    For each replicate, it runs metagen-simulator2 in parallel on each abundance table file, generating 
    metagenome sequences, and compresses the output files.

    Parameters:
    - info_folder (str): Path to the directory containing the abundance table files (each named `<sample>.tsv`).
    - out_dir (str): Path to the directory where the output for each replicate will be stored.
    - target_reads_number (int): Number of reads to simulate in each metagenome.
    - ref_genomes (str): Path to the directory containing the genomes of reference to be used.
    - read_length (int, optional): Length of each read in the simulation (default is 150).
    - replicates (int, optional): Number of replicate simulations to generate (default is 1).
    - cpu (int, optional): Number of CPU cores to use for parallel processing (default is 20).

    Dependencies:
    - `parallel`: Ensure that `parallel` is installed and accessible in the environment.
    - `pigz`: The `pigz` command should be installed for fast parallel gzip compression.
    """

    for i in range(replicates):
        # creating output folder for each replicate
        out_folder_for_replicate = os.path.join(out_dir, str(i))

        # ensuring the output directory for the replicate exists
        os.makedirs(out_folder_for_replicate, exist_ok=True)

        # building the command to find .tsv files and run metagen-simulator2 on each
        command = f"""find {shlex.quote(info_folder)} -name "*.tsv" | parallel -j {cpu} \
        metagen-simulator2 -t {read_length} -r {target_reads_number} --seed {i} --fastq {{}} {shlex.quote(ref_genomes)}"""

        # executing the simulation command
        subprocess.run(command, shell=True, check=True)

        # moving all .fastq files generated by metagen-simulator2 to the output folder for the replicate
        move_command = f"mv -v *.fastq {shlex.quote(out_folder_for_replicate)}/"
        subprocess.run(move_command, shell=True, check=True)

        # compressing the resulting .fastq files in the output folder using pigz
        compress_command = f"pigz {shlex.quote(out_folder_for_replicate)}/*.fastq"
        subprocess.run(compress_command, shell=True, check=True)

if __name__ == "__main__":
    # first generating NanoSim samples
    nanosim_simulation(info_folder="abundance_nanosim_format", out_dir="samples/nanopore",
                       nanosim_model_and_prefix="../nanosim_models/metagenome_ERR3152366_Log/training",
                       replicates=1, cpu=21, parallel=2, sample_range=(0, 31))
    # then generatin samples with metagen-simulator2
    metagen_simulation(info_folder="abundance_metagen_format", out_dir="samples/illumina",
                      ref_genomes="genomes_for_simulation",
                      target_reads_number=37000000, read_length=150, replicates=2,
                       cpu=30)