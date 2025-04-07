Here, we use species distribution as given in [*Impact of simulation and reference catalogues on the evaluation of taxonomic profiling pipelines*](https://doi.org/10.1099/mgen.0.001330) to simulate our samples. 
For three species (*E. coli*, *K. pneumoniae* and *B. infantis*, since they have the most isolate-sequenced genomes available at the UHGG) we also add
some of their strains and the total abundance of the strains equal to the original abundance of the species as given by the authors.
These strains are:

| MGnify genome | Species representative | Species                  | Note                                                                        |
| ------------- | ---------------------- | ------------------------ | --------------------------------------------------------------------------- |
| MGYG000002506 | MGYG000002506          | Escherichia coli D       | Species representative (originally present in the paper used for reference) |
| MGYG000055572 | MGYG000002506          | Escherichia coli D       | Added for the present simulations                                           |
| MGYG000177238 | MGYG000002506          | Escherichia coli D       | Added for the present simulations                                           |
| MGYG000145940 | MGYG000002506          | Escherichia coli D       | Added for the present simulations                                           |
| MGYG000262994 | MGYG000002506          | Escherichia coli D       | Added for the present simulations                                           |
| MGYG000002538 | MGYG000002538          | Klebsiella pneumoniae    | Species representative (originally present in the paper used for reference) |
| MGYG000186038 | MGYG000002538          | Klebsiella pneumoniae    | Added for the present simulations                                           |
| MGYG000163129 | MGYG000002538          | Klebsiella pneumoniae    | Added for the present simulations                                           |
| MGYG000223036 | MGYG000002538          | Klebsiella pneumoniae    | Added for the present simulations                                           |
| MGYG000131723 | MGYG000002538          | Klebsiella pneumoniae    | Added for the present simulations                                           |
| MGYG000001292 | MGYG000001292          | Bifidobacterium infantis | Species representative (originally present in the paper used for reference) |
| MGYG000147988 | MGYG000001292          | Bifidobacterium infantis | Added for the present simulations                                           |
| MGYG000066600 | MGYG000001292          | Bifidobacterium infantis | Added for the present simulations                                           |
| MGYG000006745 | MGYG000001292          | Bifidobacterium infantis | Added for the present simulations                                           |
| MGYG000074969 | MGYG000001292          | Bifidobacterium infantis | Added for the present simulations                                           |

[`build_abund_tables.py`](build_abund_tables.py) downloads the [species distrbution](https://raw.githubusercontent.com/gmtsciencedev/microbiome-pipeline-benchmarking/refs/heads/main/reference/refMet4_species_composition.tsv) given by the authors, select 25 healthy and 25 CRC-associated samples, identify the three target species,
download the MGnify genomes of all species and the four strains we want to add, recalculate abundance to adjust it for strains, produce tables for metagen-simulator2 and NanoSim.

[`simulate_samples.py`](simulate_samples.py) simulates the short and long-reads sequencing of these samples using the abundance tables.

All simulated samples were then downsized at a depth of 10 Gbp using [`downsizing_simulated_samples.py`](downsizing_simulated_samples.py).