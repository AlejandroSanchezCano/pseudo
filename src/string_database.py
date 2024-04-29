# Built-in modules
import json
import subprocess
from collections import defaultdict

# Third-party modules
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt

# Custom modules
from logger import logger
from database import Database
from interactor import Interactor

class STRING(Database):

    # Override static attributes
    database = 'STRING'
    path = Database.path + 'STRING/'
    uniprot_alias = '.aliases.v12.0.txt'

    def mads_vs_all(self) -> pd.DataFrame: # %0d doesn't work properly so we need to run multiple requests and concatenate the dataframes
        # Initialize return variable
        string_df = pd.DataFrame()

        # Populate dataframe with individual requests
        geneList = [string_id for interactor in self.in_database() for string_id in interactor.string_id]
        for string_id in tqdm(geneList):
            # STRING network
            cmd = f'curl "https://string-db.org/api/json/network?identifiers={string_id}&required_score=0"'
            string_results = subprocess.run(cmd, capture_output = True, text = True, shell = True).stdout
            string_json = json.loads(string_results)
            # Concatenate with STRING dataframe network
            df = pd.DataFrame.from_dict(string_json, orient='columns')
            string_df = pd.concat([string_df, df], ignore_index = True)
        
        # Remove duplicated rows
        string_df = string_df.drop_duplicates(ignore_index = True)

        # Logging
        logger.info(f'{len(string_df)} MADS-vs-all interactions retrieved from {len(geneList)} UniProt IDs in {self.database}')

        return string_df
    
    def mads_vs_mads(self) -> pd.DataFrame: 

        # MADS-vs-all
        string_df = self.mads_vs_all()

        # Filter MADS-vs-MADS
        interactor_list = [interactor for interactor in self.in_database()]
        mads_string_ids = set([string_id for interactor in interactor_list for string_id in interactor.string_id])
        mads_interactions = string_df.apply(lambda df: set([df.stringId_A, df.stringId_B]) - mads_string_ids == set(), axis = 1)
        string_df = string_df[mads_interactions]

        # Logging
        logger.info(f'{len(string_df)} MADS-vs-MADS interactions retrieved from {len(interactor_list)} UniProt IDs in {self.database}')

        return string_df
    
    def filter_score(self, string_df: pd.DataFrame, **kwargs: dict[str, int]) -> pd.DataFrame:
        for score, value in kwargs.items():
            string_df = string_df.loc[string_df[score] >= value]

        return string_df
    
    def size_per_filter(self, df: pd.DataFrame) -> None:
        scores = [column for column in df.columns if 'score' in column]
        cutoffs = np.arange(0, 1.01, 0.01)

        score_to_size = defaultdict(list)
        for score in tqdm(scores):
            for cutoff in cutoffs:
                network = self.filter_score(df, **{score: cutoff})
                score_to_size[score] += [len(network)]
        score_to_size = pd.DataFrame(score_to_size, index = cutoffs)
        sns.lineplot(data = score_to_size, dashes= False)
        plt.xlabel('Score cutoff')
        plt.ylabel('Number of interactions')
        plt.show()

    def scores_in_taxa(self, df: pd.DataFrame, taxa = '3702') -> None:
        scores = [column for column in df.columns if 'score' in column]
        df = df[df['ncbiTaxonId'] == '3702'][scores]
        sns.violinplot(data = df)
        plt.xlabel('STRING scores')
        plt.ylabel(f'Score value')
        plt.show()

    def score_per_taxa(self, df: pd.DataFrame, score = 'score') -> None:
        order = df.groupby('ncbiTaxonId')[score].median().sort_values(ascending = False).index
        sns.violinplot(data = df, x = 'ncbiTaxonId', y = score, order = order)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.xlabel('Species')
        plt.ylabel(f'STRING {score}')
        plt.show()

    def translate_ids(self, df: pd.DataFrame) -> pd.DataFrame:

        mads_uniprot_ids = [interactor.uniprot_id for interactor in self._iterate_interactors()]
        _, string_to_uniprot = self._get_aliases()

        def to_uniprot(string_id):

            uniprot_ids = string_to_uniprot[str(string_id)]

            if len(uniprot_ids) == 1:
                return uniprot_ids[0]
            
            merged_uniprot_ids = [uniprot_id for uniprot_id in uniprot_ids if uniprot_id in mads_uniprot_ids]

            if len(merged_uniprot_ids) == 1:
                return merged_uniprot_ids[0]
            
            swissprot_uniprot_ids = [uniprot_id for uniprot_id in merged_uniprot_ids if Interactor.unpickle(uniprot_id).uniprot_info['Section'] == 'Swiss-Prot']

            if len(swissprot_uniprot_ids) == 1:
                return swissprot_uniprot_ids[0]
            
            logger.error(f'STRING ID {string_id} cannot be translated properly')
        
        df['UNIPROT_ID_A'] = df['stringId_A'].apply(to_uniprot)
        df['UNIPROT_ID_B'] = df['stringId_B'].apply(to_uniprot)   

        return df


    
            
if __name__ == '__main__':
    s = STRING()
    #s.map_ids()
    #s.mads_vs_mads()
    df = pd.read_pickle('s.pkl')

    # Analyze scores
    #s.size_per_filter(df)
    #s.scores_in_taxa(df, taxa = '3702')
    #s.score_per_taxa(df, score = 'score')
    #s.score_per_taxa(df, score = 'escore')

    #df = s.translate_ids(df)

'''
import utils
a = utils.unpickling('/home/asanche2/asanche2/ppi/data/Networks/STRING_MADS vs MADS.pkl')
s.size_per_filter(a)
'''
