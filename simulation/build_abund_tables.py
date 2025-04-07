import pandas as pd
import sys
import os
import subprocess
import json
import numpy as np
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import get_uhgg_genomes_list as gugl

class forSimulation:
    def __init__(self, samples_metadata: list):
        # loading abundance data
        print("Loading article's abundances data...")
        self.abundance_data = self._load_abundance_data("https://raw.githubusercontent.com/gmtsciencedev/microbiome-pipeline-benchmarking/refs/heads/main/reference/refMet4_species_composition.tsv")
        
        # loading UHGG metadata
        print("Loading UHGG gut metadata...")
        self.uhgg_metadata = self._load_uhgg_metadata()
        
        # merging the metadata with abundance data
        print("Merging data...")
        self.data = self._merge_metadata_abundance(self.abundance_data, self.uhgg_metadata)
        
        # extracting the sample names (columns of the abundance data)
        self.samples = list(self.abundance_data.columns)

        # loading samples metadata
        self.samples_metadata = self._load_samples_metadata(samples_metadata)

    @staticmethod
    def _load_abundance_data(path: str) -> pd.DataFrame:
        """
        Loads the abundance data from a TSV file.
        
        Parameters:
        - path (str): Path to the TSV file with abundance data.
        
        Returns:
        - pd.DataFrame: DataFrame with species abundance data.
        """
        return pd.read_csv(path, sep='\t', index_col=0)

    @staticmethod
    def _load_uhgg_metadata() -> pd.DataFrame:
        """
        Loads the UHGG metadata.
        
        Returns:
        - pd.DataFrame: DataFrame with UHGG metadata indexed by genome.
        """
        metadata = pd.read_csv("https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv",
                                sep='\t')
        
        return metadata.set_index('Genome')

    @staticmethod
    def _merge_metadata_abundance(abundances: pd.DataFrame, uhgg_metadata: pd.DataFrame) -> pd.DataFrame:
        """
        Merges UHGG metadata with species abundance data on genome index.
        
        Parameters:
        - abundances (pd.DataFrame): DataFrame with species abundance data.
        - uhgg_metadata (pd.DataFrame): DataFrame with UHGG metadata.
        
        Returns:
        - pd.DataFrame: Merged DataFrame with UHGG metadata and abundance data.
        """
        return pd.merge(uhgg_metadata[['Lineage']], abundances, left_index=True, right_index=True, how='right')

    @staticmethod
    def _divide_values(values: list, n: int):
        """
        Divides each value in the list 'values' into 'n' unequal parts simultaneously.

        Parameters:
        - values (List[float]): A list of values to be divided.
        - n (int): The number of parts to divide each value into.

        Returns:
        - List[List[float]]: A list of lists where each inner list contains n parts that sum to the corresponding value in 'values'.
        """
        values = np.array(values)
        
        # generating random proportions for each value in the shape (len(values), n)
        random_proportions = np.random.random((len(values), n))
        
        # normalizing each row to sum to 1
        normalized_proportions = random_proportions / random_proportions.sum(axis=1, keepdims=True)
        
        # getting the new abundances
        parts = normalized_proportions * values[:, np.newaxis]
        
        return parts.tolist()

    def _load_samples_metadata(self, path_to_metadata: list) -> pd.DataFrame:
        """ 
        Loads and merges multiple sample metadata files, filters for unique entries, 
        and keeps only samples present in the author's samples list.
        
        Parameters:
        - path_to_metadata (list): List of paths to TSV files, each indexed by sample with a "phenotype" column.
        
        Returns:
        - pd.DataFrame: Merged metadata DataFrame with unique samples, indexed by sample.
        """
        
        list_metadata = []

        for path in path_to_metadata:
            list_metadata.append(
                pd.read_csv(path, sep='\t', index_col=0)
            )

        metadata = pd.concat(list_metadata)
        metadata = metadata.reset_index()
        metadata = metadata.drop_duplicates()
        metadata = metadata.set_index('biosample_accession')

        # keeping only samples that appear in authors table
        return metadata.loc[self.samples]
    
    def get_list_classes(self):
        """
        Retrieves a list of unique classes (phenotypes) from the samples metadata.
        
        Returns:
        - List[str]: A list of unique phenotypes present in the samples metadata.
        """
        return self.samples_metadata['phenotype'].unique().tolist()
    
    def select_samples(self, class_a: str, class_b: str, samples_by_class: int, seed: int = 999):
        """
        Selects a specified number of samples from two given classes (phenotypes) 
        and returns their indices sorted in ascending order.

        Parameters:
        - class_a (str): The first phenotype class from which to select samples.
        - class_b (str): The second phenotype class from which to select samples.
        - samples_by_class (int): The number of samples to select from each class.
        - seed (int, optional): Random seed for reproducibility.

        Returns:
        - list: A sorted list of selected samples from both classes.
        """
        # selecting samples of class_a
        search_criterion_a = self.samples_metadata['phenotype'] == class_a
        samples_a = self.samples_metadata[search_criterion_a].sample(n=samples_by_class, random_state=seed).index

        # selecting samples of class_b
        search_criterion_b = self.samples_metadata['phenotype'] == class_b
        samples_b = self.samples_metadata[search_criterion_b].sample(n=samples_by_class, random_state=seed).index

        return list(np.sort(np.concatenate((samples_a, samples_b))))
    
    def get_potential_strains(self, n: int = 10):
        """ 
        Returns a DataFrame with UHGG genomes potentially interesting species because they contain a lot of isolate-sequenced strains

        It only uses the species given in the paper and identifies the ones having the more isolate genomes in the
        UHGG. It outputs top-n species
        """

        species_in_simulation = list(self.data.index)
        uhgg_metadata_filtered_on_species = self.uhgg_metadata[self.uhgg_metadata['Species_rep'].isin(species_in_simulation)]
        # only keeping isolate genomes
        uhgg_metadata_filtered_on_species = uhgg_metadata_filtered_on_species[uhgg_metadata_filtered_on_species['Genome_type'] == 'Isolate']
        # get number of genomes by species to identify interesting species to use
        genomes_by_species = uhgg_metadata_filtered_on_species['Species_rep'].value_counts().head(n).to_frame()

        return pd.merge(genomes_by_species, self.uhgg_metadata[['Lineage']], how='left', left_index=True, right_index=True)
    
    def sample_strains_genomes(self, species: str, n: int, seed: int = 999):
        """ 
        Randomly select genomes of isolates of a given species.
        
        This method filters the UHGG metadata to retrieve strains belonging to 
        the specified species and randomly samples a specified number of 
        genomes, excluding the representative genome of that species.
        
        Parameters:
        - species (str): The species for which genomes of isolates will be sampled.
        - n (int): The number of genomes to randomly select.
        - seed (int, optional): Random seed for reproducibility.
        
        Returns:
        - List[str]: A list of UHGG genomes id.
        """
        
        uhgg_metadata_filtered_on_strains = self.uhgg_metadata[(self.uhgg_metadata['Species_rep'] == species) & (self.uhgg_metadata['Genome_type'] == 'Isolate')]
        # don't include the genome of the species representative itself since it is already included in the abundances
        uhgg_metadata_filtered_on_strains = uhgg_metadata_filtered_on_strains.drop(index=species)

        return list(uhgg_metadata_filtered_on_strains.sample(n, random_state=seed).index)
    
    def recalculate_abundances_using_strains(self, strains: dict, seed: int = 999):
        """ 
        Recalculates abundances by dividing the abundance of specific species 
        across strains within each sample.

        Parameters:
        - strains (dict): A dictionary where each key is a species representative and each value 
                          is a list of genome identifiers representing the strains for that species.
                          For example, {'B': ['B1', 'B2']} means species 'B' has two strains, 'B1' and 'B2' (+ 'B'
                          the original genome in the article).
        - seed (int): Random seed for reproducibility in abundance division.

        Process:
        - For each species that has associated strains in the `strains` dictionary, 
          this method recalculates the abundance values by dividing the species' original abundance 
          across its strains and the original species representative.
        - The sum of the new abundances for the species and its strains equals the original abundance.

        Example:
        Given an initial abundance data like this:
        
        | species | sample1 | sample2 |
        |---------|---------|---------|
        | A       | 0.4     | 0       |
        | B       | 0.2     | 0.3     |
        | C       | 0       | 0.1     |
        | D       | 0.4     | 0.7     |

        And a strains dictionary: {'B': ['B1', 'B2']}
        
        This method might yield:
        
        | species | sample1    | sample2    |
        |---------|------------|------------|
        | A       | 0.4        | 0          |
        | B       | 0.05366598 | 0.115501   |
        | B1      | 0.07232264 | 0.03517515 |
        | B2      | 0.07401139 | 0.04932385 |
        | C       | 0          | 0.1        |
        | D       | 0.4        | 0.7        |

        Note:
        The new abundances satisfy the condition:
        `B (original) = B (new) + B1 + B2` for all samples.

        Returns:
        - pd.DataFrame: A DataFrame with updated abundances for each species and its strains.
        """
        
        # copying the abundance data to avoid modifying the original dataset
        abundances_data = self.abundance_data.copy()

        for species in strains:
            # getting the abundance data for the species across all samples
            abundance_in_sample = abundances_data.loc[species].to_list()
            
            # diving abundance across the strains and the original species representative
            # these new abundances are not equal between strains in a same sample
            new_abundances = self._divide_values(
                abundance_in_sample, 
                len(strains[species]) + 1  # +1 for the original species representative
            )
            
            # transposing and preparing abundance values for each genome (original + strains)
            new_abundances = list(zip(*new_abundances))
            list_genomes = [species] + strains[species]  # including the species and its strains in the index
            
            # creating a DataFrame with recalculated abundances for each strain and species representative
            new_abundances_df = pd.DataFrame(new_abundances, columns=self.samples, index=list_genomes)

            # replacing the original species row with new abundances for this species and strains
            abundances_data = abundances_data.drop(index=species)  # removing the old row for the species
            abundances_data = pd.concat([abundances_data, new_abundances_df])  # adding recalculated rows

        return abundances_data
    
    def download_and_process_genomes_from_uhgg(self, abundance_data: pd.DataFrame, download_dir: str, sleep: int = 2):
        """ 
        Downloads and processes genomes used in these simulations.

        This method merges genome download links with the abundance data, downloads the relevant genome files, 
        and processes them by extracting genomic sequences from .GFF files.

        Parameters:
        - abundance_data (pd.DataFrame): A DataFrame containing abundance data for each genome and samples.
        - download_dir (str): Path to the directory where genomes will be downloaded.
        - sleep (int): Delay (in seconds) between downloads.
        """

        # adding FTP download link to abundance data
        abundance_data = pd.merge(abundance_data, self.uhgg_metadata[['FTP_download']], 
                                  left_index=True, right_index=True, how='left')
        
        # downloading the genomes
        gugl.download_genomes(abundance_data, download_dir, sleep)

        # processing the genomes (extracting sequences from .GFF)
        gugl.processing_uhgg_genomes(download_dir)

    @staticmethod
    def build_metagen_simulator2_tables(abundance_data: pd.DataFrame, output_dir: str):
        """
        Converts abundance data into tables compatible with the metagen-simulator2 program.

        This method takes a DataFrame of abundance data, saves it temporarily as a TSV file, and then runs an external 
        script to convert this file into the format required by the metagen-simulator2 tool. The converted tables 
        are saved in the specified output directory.

        Parameters:
        - abundance_data (pd.DataFrame): DataFrame containing abundance data where each row represents a genome 
          and each column represents a sample.
        - output_dir (str): Path to the directory where the output tables for metagen-simulator2 will be saved.

        Dependencies:
        - Requires a Python script (`build_metagen_simulator_tables.py`) available at the specified path, which handles 
          the conversion of the TSV file into metagen-simulator2-compatible tables.
        """

        # exporting the abundance data table for temporary use
        abundance_data.index.name = 'Genome'
        abundance_data.to_csv("tmp.tsv", sep='\t')

        conversion_script = "../build_metagen_simulator_tables.py"
        cmd = f"python3 {conversion_script} --out_dir {output_dir} tmp.tsv"

        # running the `build_metagen_simulator_tables.py` command in the shell
        subprocess.run(cmd, shell=True)

        # removing the temporary file after the conversion process
        os.remove("tmp.tsv")

    @staticmethod
    def build_nanosim_tables(abundance_data: pd.DataFrame, output_dir: str, genomes_folder: str, target_number_reads: int):
        """ 
        Converts abundance data into tables compatible with the NanoSim program. 

        This method exports abundance data to a temporary TSV file and then runs an external script to convert 
        the TSV file into the format required by the NanoSim tool. The number of reads to be simulated afterwards must 
        be specified as `target_number_reads`.

        Parameters:
        - abundance_data (pd.DataFrame): DataFrame containing abundance data, where each row represents a genome 
          and each column represents a sample.
        - output_dir (str): Path to the directory where the output tables for NanoSim will be saved.
        - genomes_folder (str): Path to the folder containing the genome files to be used by NanoSim.
        - target_number_reads (int): The number of reads to simulate, as required by the NanoSim program.

        Dependencies:
        - Requires a Python script (`build_nanosim_tables.py`) available at the specified path, which handles 
          the conversion of the TSV file into NanoSim-compatible tables.
        """

        # exporting the abundance data table for temporary use
        abundance_data.index.name = 'Genome'
        abundance_data.to_csv("tmp.tsv", sep='\t')

        conversion_script = "../build_nanosim_tables.py"
        cmd = f"python3 {conversion_script} --reads {int(target_number_reads)} --out_dir {output_dir} tmp.tsv {genomes_folder}"

        # running the `build_metagen_simulator_tables.py` command in the shell
        subprocess.run(cmd, shell=True)

        # removing the temporary file after the conversion process
        os.remove("tmp.tsv")

if __name__ == '__main__':
    sim_data = forSimulation(samples_metadata=["PRJEB6070.tsv", "PRJEB7774.tsv"])

    # interesting species we could use to introduce strains in the simulation
    species_with_lot_of_strains = sim_data.get_potential_strains()

    #                count                                                                                                                                            Lineage
    # Species_rep                                                                                                                                                            
    # MGYG000002506   3922                d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli_D
    # MGYG000002538    344              d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Klebsiella;s__Klebsiella pneumoniae
    # MGYG000002369    258           d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Clostridioides;s__Clostridioides difficile
    # MGYG000002353    214                             d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus_B;s__Enterococcus_B faecium
    # MGYG000001292    146           d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium infantis
    # MGYG000002478    143                                  d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Phocaeicola;s__Phocaeicola dorei
    # MGYG000001346    112                              d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides uniformis
    # MGYG000003683     70  d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium pseudocatenulatum
    # MGYG000002438     69                     d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Tannerellaceae;g__Parabacteroides;s__Parabacteroides distasonis
    # MGYG000002395     68       d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Actinomycetales;f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium adolescentis

    # selecting 4 more Escherichia coli D strains
    e_coli_strains = sim_data.sample_strains_genomes('MGYG000002506', 4)

    # selecting 4 more Klebsiella pneumoniae strains
    k_pneumoniae_strains = sim_data.sample_strains_genomes('MGYG000002538', 4)

    # not selecting 4 more Clostridioides difficile strains
    # there are a lot of C. difficile strains but it is only present in 0,87% of the samples...
    # c_difficile_strains = sim_data.sample_strains_genomes('MGYG000002369', 4)

    # selecting 4 more Bifidobacterium infantis strains
    b_infantis_strains = sim_data.sample_strains_genomes('MGYG000001292', 4)

    # constructing a dictionary of species -> strains to recalculate the abundances
    species_strains = {
        'MGYG000002506': e_coli_strains,
        'MGYG000002538': k_pneumoniae_strains,
        'MGYG000001292': b_infantis_strains
    }

    # saving the genomes id into a file
    file_name = 'strains_genomes.json'
    with open(file_name, 'w') as json_file:
        json.dump(species_strains, json_file, indent=4)

    # randomly selecting 25 healthy and 25 sick patient' samples
    selecting_samples = sim_data.select_samples('D006262', # health
                                                'D015179', # neoplasms
                                                samples_by_class=25)

    # introducing these strains in the original abundances given in the article
    # E. coli, K. pneumoniae and C. difficile abundances are still the same in samples but they now have strains
    # (= strains abundances sum to the original species abundance)
    new_abundances = sim_data.recalculate_abundances_using_strains(species_strains)
    
    # from now we will only use the random samples we got
    new_abundances = new_abundances[selecting_samples]

    # downlading the genomes from the UHGG server
    sim_data.download_and_process_genomes_from_uhgg(new_abundances, 'genomes_for_simulation')

    # converting the abundances into tables for use by metagen-simulator2
    sim_data.build_metagen_simulator2_tables(new_abundances, 'abundance_metagen_format')

    # converting the abundances into tables for use by NanoSim
    sim_data.build_nanosim_tables(new_abundances, 'abundance_nanosim_format', 'genomes_for_simulation',
                                  target_number_reads=2.7e6)
