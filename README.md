## Comparative Evaluation of Short-Read, Long-Read, and Hybrid Assemblies for MAG Recovery in Human Fecal Metagenomes

### Pipeline 

[Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline in [`metageno-pipeline`](metageno-pipeline/) has been used for generating most data of the study.
It has been tested with Snakemake v8.24.

The configuration file [`config/simulation.yaml`](config/simulation.yaml) contains the set of assemblers, binners and output params we used for the study.

You will need a TSV with the path to reads, etc.
Check [`metageno-pipeline/data/config_data.tsv`](metageno-pipeline/data/config_data.tsv).

To run the pipeline:

```sh
cd metageno-pipeline/

snakemake --use-conda --jobs <n> \
    --rerun-triggers mtime \
    --keep-going \
    --configfile /path/to/config/simulation.yaml
```

We just need the following output in the [`Snakefile`](metageno-pipeline/workflow/Snakefile):

```python
rule all:
        input:
            # assembly part
            expand("results/03_assembly/{assembler}/{sample}/assembly.fa.gz", 
                   assembler=ASSEMBLER + HYBRID_ASSEMBLER, sample=SAMPLES),
            expand("results/03_assembly/{assembler_long_read}/{sample_lr}/assembly.fa.gz", 
                   assembler_long_read=LONG_READ_ASSEMBLER, sample_lr=SAMPLES_LR),
            # assembly qc
            expand("results/04_assembly_qc/quast/{assembler}/{sample}/combined_reference/report.tsv", 
                   assembler=ASSEMBLER + HYBRID_ASSEMBLER, sample=SAMPLES),
            expand("results/04_assembly_qc/quast/{assembler_lr}/{sample_lr}/combined_reference/report.tsv", 
                   assembler_lr=LONG_READ_ASSEMBLER, sample_lr=SAMPLES_LR),
            # non redundant gene catalog
            expand("results/04_assembly_qc/gene_clustering/{assembler}/non_redundant_gene_catalog.fna.gz",
                   assembler = ASSEMBLER + HYBRID_ASSEMBLER + LONG_READ_ASSEMBLER),
            # binning (step 05) + bins qc (step 06)
            expand("results/05_binning/{binner}/bins/{assembler}/{sample}",
                   binner = SHORT_READ_BINNER, assembler=ASSEMBLER + HYBRID_ASSEMBLER, sample=SAMPLES),
            expand("results/05_binning/{binner_lr}/bins/{assembler_lr}/{sample_lr}",
                   binner_lr = LONG_READ_BINNER, assembler_lr=ASSEMBLER + HYBRID_ASSEMBLER, sample_lr=SAMPLES),
            expand("results/06_binning_qc/checkm2/samples/{sample}/all_quality_reports.pdf",
                   sample=SAMPLES),
            # bins post-processing
            expand("results/08_bins_postprocessing/dRep/{ani}/{assembler}",
                   assembler=ASSEMBLER + LONG_READ_ASSEMBLER + HYBRID_ASSEMBLER, ani = ANI_THRESHOLD),
            expand("results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/{ani}/{assembler}/bins",
                   assembler=ASSEMBLER + HYBRID_ASSEMBLER + LONG_READ_ASSEMBLER, ani = ANI_THRESHOLD),
```

For Skani pairwise comparison of produced MAGs use the supplementary script [`metageno-pipeline/workflow/scripts/other_scripts/skani_analysis.py`](metageno-pipeline/workflow/scripts/other_scripts/skani_analysis.py).

If you want to have assembly metrics based on only specific reference genomes introduced in each simulated samples (NA50, NGA50... only on genomes of simulation) better use the file [`scripts/metaquast_commands.txt`](scripts/metaquast_commands.txt) that already contains the exact reference genomes by sample for MetaQUAST.

### Scripts

ANI comparison to references before and after dereplication script is at [`scripts/skani_bins_vs_references.sh`](scripts/skani_bins_vs_references.sh).
It will save the results in `scripts/skani_search_against_ref`.

### Data

#### Simulated dataset

Code for simulation is stored in [`simulation`](simulation/).
We've used metagen-simulator2 for short-read simulation and Nanosim v3.2.2 for long-read simulation.

[`simulation/README.md`](simulation/README.md) provides more information.

#### Results

Excel table with stats on simulated samples, assembly and binning is available at [`data/supplementary_materials.xlsx`](data/supplementary_materials.xlsx).

| Sheet                         | Description                                                                                                                                                                                                                                                                                                                              |
| :---------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Simulated samples             | Abundances used by sample for the simulation                                                                                                                                                                                                                                                                                             |
| Sequencing and Assembly Stats | Information on simulated samples (total length, number of reads), metrics on assembly (number of contigs, assembly size, min/max lenght of contigs, N50, L50, NA50, LA50, auNA, $\frac{\text{N50} - \text{NA50}}{\text{N50}}$, $\frac{\text{L50} - \text{LA50}}{\text{L50}}$), samples estimated richness by Meteor (gene and MSP count) |
| Genomes NGA50 LGA50           | Value of genomes' NGA50 and LGA50 by assembly approach                                                                                                                                                                                                                                                                                   |
| Bins Stats Before dRep        | CheckM2 metrics of SemiBin2 bins (before dereplication)                                                                                                                                                                                                                                                                                  |
| MAGs Stats (dRep 97% ANI)     | CheckM2 metrics of SemiBin2 bins (after dereplication and filtration on completeness and contamination step)                                                                                                                                                                                                                             |
| Bins Before dRep (Skani)      | Skani search of bins against genomes of reference (before dereplication)                                                                                                                                                                                                                                                                 |
| = (confusion matrix)      | Confusion matrix of skani search of bins against genomes of reference, with accuracy and F1 score calculated (before dereplication)                                                                                                                                                                                                                                                                 |
| Bins After dRep (Skani)       | Skani search of bins against genomes of reference (after dereplication and filtration on completeness and contamination step)                                                                                                                                                                                                            |

These tables are also available as individual TSV in [`data/tables`](data/tables).