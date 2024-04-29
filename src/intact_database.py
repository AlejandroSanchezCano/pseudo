# Built-in modules
import subprocess
from io import StringIO
from typing import Literal
from collections import defaultdict

# Third-party modules
import pandas as pd
from tqdm import tqdm

# Custom modules
from logger import logger
from database import Database

# TODO: we might need to also parse the negative file

class IntAct(Database):

    # Override static attributes
    database = 'IntAct'
    path = Database.path + 'IntAct/'
    uniprot_alias = 'links.dat'

    @staticmethod
    def _get_aliases() -> tuple[dict[str, list[float]], None]:
        '''
        IntAct uses UniProt IDs as main IDs, so despite parsing a file that maps
        the IDs to the number of interactions in IntAct, we only need the list 
        of IDs with information in the database.

        Returns
        -------
        dict[str, list[float]]
            Maps UniProt ID to number of interactions.
        '''

        # Initialize return variable
        aliases_dict = defaultdict(int)

        # Open database file
        database_file = f'{IntAct.path}uniprot{IntAct.uniprot_alias}'
        with open(database_file, 'r') as handle:

            # Parse database file
            for line in tqdm(handle.read().splitlines()):
                _, uniprot_id, n_interactions = line.split('; ')

                # Add to alias-mapping dictionary
                aliases_dict[uniprot_id] = float(n_interactions)
        
        # Logging
        logger.info(f'{len(aliases_dict)} total UniProt IDs mapped to {IntAct.database}')

        return aliases_dict, None
    
    def mads_vs_all(self, type: Literal['+', '-'] = '+') -> pd.DataFrame:
        '''
        Parses the IntAct PPI file for positive or negative interactions and
        retrieves the MADS vs. all PPIs. 

        Returns
        -------
        pd.DataFrame
            MADS vs. all data frame network.
        '''
        # Find MADS UniProt IDs in IntAct file
        mads_interactors_ids = [interactor.intact_id for interactor in self.in_database()]
        type_network = '' if type == '+' else '_negatives'
        file_path = f'/home/asanche2/asanche2/ppi/data/IntAct/intact{type_network}.txt'
        cmd = f'cat {file_path} | grep -iE "^#ID|{"|".join(mads_interactors_ids)}"'
        intact_results = subprocess.run(cmd, capture_output = True, text = True, shell = True).stdout

        # No MADS in - PPIs
        if intact_results:
            df = pd.read_csv(StringIO(intact_results), sep = '\t')
        else:
            df = pd.DataFrame()

        # Logging
        logger.info(f'{len(df)} {type} MADS-vs-all rows retrieved from {len(mads_interactors_ids)} UniProt IDs in {self.database}')

        return df
    
    def mads_vs_mads(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Filters the MADS-vs-all PPI network to contain only MADS-vs-MADS PPIs.
        There are Interactors whose main ID is not UniProt, but only one is a
        MADS protein:
        - EBI-15199253 -> P17839 (AG)
        - At1g28710 NO MADS
        - At3g26140 NO MADS
        - At3g21090 NO MADS
        - At1g06080 NO MADS
        - EBI-15192079 NO MADS

        Returns
        -------
        pd.DataFrame
            MADS-vs-MADS interaction network.
        '''
        # Change EBI-15199253 to AG (P17839)
        df.loc[df['#ID(s) interactor A'] == 'intact:EBI-15199253', ['#ID(s) interactor A']] = 'uniprotkb:P17839'
        df.loc[df['ID(s) interactor B'] == 'intact:EBI-15199253', ['ID(s) interactor B']] = 'uniprotkb:P17839'

        # Filter MADS-vs-MADS
        interactor_list = [interactor for interactor in self.in_database()]
        mads_intact_ids = set([interactor.intact_id for interactor in interactor_list])
        df['A'] = df['#ID(s) interactor A'].apply(lambda x: x.split(':')[1])
        df['B'] = df['ID(s) interactor B'].apply(lambda x: x.split(':')[1])
        mads_interactions = df.apply(lambda df: set([df.A, df.B]) - mads_intact_ids == set(), axis = 1)
        df = df[mads_interactions]

        # Logging
        logger.info(f'{len(df)} MADS-vs-MADS rows retrieved from {len(interactor_list)} UniProt IDs in {self.database}')

        return df
    
    def standarize(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Format data frame interaction network to accommodate standard naming
        convention of columns to homogenize the data frames from different 
        databases.

        Parameters
        ----------
        df : pd.DataFrame
            Data frame interaction network.

        Returns
        -------
        pd.DataFrame
            Standarize data frame interaction network.
        '''
        df['A-B'] = df.apply(lambda x: '-'.join(sorted([x['A'], x['B']])), axis=1)

        return df

if __name__ == '__main__':
    '''Test class'''
    intact = IntAct()
    intact.map_ids()
    df = intact.mads_vs_all('-')
    df = intact.mads_vs_all('+')
    df = intact.mads_vs_mads(df)
    df = intact.standarize(df)
    intact.save(df, 'MADS')


