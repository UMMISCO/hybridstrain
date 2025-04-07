#!/bin/bash

module purge
# source activate <conda_env_with_skani>/skani.0.2.2

# building Skani database for reference genomes (given genomes are in metageno-pipeline/data/reference_genomes)
skani sketch metageno-pipeline/data/reference_genomes/*.fa * -o database

# before dereplication (bins)
# metaflye
skani search metageno-pipeline/results/05_binning/semibin2/bins/metaflye/*/bins/*.fa.gz \
    -d database/ -o "scripts/skani_search_against_ref/before_dereplication/metaflye_res.tsv" \
    -t 20
 
# megahit
skani search metageno-pipeline/results/05_binning/semibin2/bins/megahit/*/bins/*.fa.gz \
    -d database/ -o "scripts/skani_search_against_ref/before_dereplication/megahit_res.tsv" \
    -t 20
 
# hybridspades
skani search metageno-pipeline/results/05_binning/semibin2/bins/hybridspades/*/bins/*.fa.gz \
    -d database/ -o "scripts/skani_search_against_ref/before_dereplication/hybridspades_res.tsv" \
    -t 20

# after dereplication (MAGs)
# metaflye
skani search metageno-pipeline/results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/97/metaflye/bins/*.fa \
    -d database/ -o "scripts/skani_search_against_ref/metaflye_res.tsv" \
    -t 20
 
# megahit
skani search metageno-pipeline/results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/97/megahit/bins/*.fa \
    -d database/ -o "scripts/skani_search_against_ref/megahit_res.tsv" \
    -t 20
 
# hybridspades
skani search metageno-pipeline/results/08_bins_postprocessing/dereplicated_genomes_filtered_by_quality/97/hybridspades/bins/*.fa \
    -d database/ -o "scripts/skani_search_against_ref/hybridspades_res.tsv" \
    -t 20