import os

def main():
    # 1: listing all files in "abundance_metagen_format/" and extracting their name, without the extension
    folder = "abundance_metagen_format/"
    sample_ids = [os.path.splitext(f)[0] for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]

    # listing to store all commands for GNU Parallel
    commands = []

    # 2: generating the commands for each simulated sample
    for sample in sample_ids:
        for i in range(3):
            # ensure output directories exist
            os.makedirs(f"samples_downsized/nanopore/{i}", exist_ok=True)
            os.makedirs(f"samples_downsized/illumina/{i}", exist_ok=True)

            # Nanopore command
            nanopore_cmd = (
                f"rasusa reads --bases 10g "
                f"-o samples_downsized/nanopore/{i}/{sample}_nanopore.fasta.gz "
                f"samples/nanopore/{i}/{sample}_nanopore.fasta.gz"
            )
            commands.append(nanopore_cmd)

            # Illumina command
            illumina_cmd = (
                f"rasusa reads --bases 10g "
                f"-o samples_downsized/illumina/{i}/{sample}_1.fastq.gz "
                f"-o samples_downsized/illumina/{i}/{sample}_2.fastq.gz "
                f"samples/illumina/{i}/{sample}_1.fastq.gz "
                f"samples/illumina/{i}/{sample}_2.fastq.gz"
            )
            commands.append(illumina_cmd)

    # 3: writing the commands to a file for GNU Parallel
    with open("commands.txt", "w") as f:
        f.write("\n".join(commands))

    # 4: running the commands using GNU Parallel with a parallelism of 45
    parallel_cmd = "parallel -j 54 < commands.txt"
    print(f"Executing all commands with: {parallel_cmd}")
    os.system(parallel_cmd)

if __name__ == "__main__":
    main()
