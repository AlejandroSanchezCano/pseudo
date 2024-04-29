# Built-in modules
import subprocess
from io import StringIO

# Third-party modules
import pandas as pd
from tqdm import tqdm
from biogridpy.biogrid_client import BioGRID

# Custom modules
from logger import logger
from database import Database
from interactor import Interactor

# TODO: STILL BIOGRID NEEDS A WAY TO UNEXPLODE THE methods part -> groupby('A-B').agg(list)

class BioGrid(Database):

    # Override static attributes
    database = 'BioGRID'
    path = Database.path + 'BioGRID/'
    uniprot_alias = '-IDENTIFIERS-4.4.230.tab.txt'

    @classmethod
    def remove_empty_accessions(cls) -> None:
        '''
        It is not rare in BioGRID to find empty accessions, where a protein has
        an associated BioGRID ID but no PPIs. Given that they carry no
        information, they will are removed at the Interactor level by using the
        BioGRID API to retrieve the interactions of a given UniProt ID and check
        whether the response is Falsy. 
        I takes 5 min for 657 accessions.
        '''
        # Count empty accessions
        counter = 0

        for interactor in tqdm(cls.in_database()):
            # Access BioGRID API
            cmd = f'curl "https://webservice.thebiogrid.org/interactions/?geneList={interactor.uniprot_id}&accesskey=b77211bccdd281cf14cca4e39926efa9"'
            response = subprocess.run(cmd, capture_output = True, text = True, shell = True).stdout

            # Disregard empty accessions
            if not response:
                counter += 1
                interactor.biogrid_id = ''
                interactor.pickle()
                
        # Logging
        logger.info(f'{counter} MADS UniProt IDs without interactions in BioGRID')

    def mads_vs_all(self) -> pd.DataFrame:
        '''
        Uses the BioGRID API wrapper to retrieve the MADS-vs-all PPI network.

        Returns
        -------
        pd.DataFrame
            MADS-vs-all interaction network.
        '''
        # Network
        BG = BioGRID(config_filepath='/home/asanche2/asanche2/miniconda3/envs/mads/lib/python3.10/site-packages/biogridpy/biogridpyrc')
        geneList = [interactor.uniprot_id for interactor in self.in_database()]
        bg_results = BG.interactions('json', geneList = geneList)
        bg_json = StringIO(bg_results.toDataFrame())
        biogrid_df = pd.read_json(bg_json, orient='index')

        # Logging
        logger.info(f'{len(biogrid_df)} MADS-vs-all rows retrieved from {len(geneList)} UniProt IDs in {self.database}')

        return biogrid_df
    
    def mads_vs_mads(self) -> pd.DataFrame: 
df['INTERACTION_ID'] = df.apply(lambda x: '-'.join(sorted([x['A'], x['B']])), axis=1)
        # MADS-vs-all
        biogrid_df = self.mads_vs_all()

        # Filter MADS-vs-MADS
        interactor_list = [interactor for interactor in self.in_database()]
        mads_biogrid_ids = set([int(biogrid_id) for interactor in interactor_list for biogrid_id in interactor.biogrid_id])
        mads_interactions = biogrid_df.apply(lambda df: set([df.BIOGRID_ID_A, df.BIOGRID_ID_B]) - mads_biogrid_ids == set(), axis = 1)
        biogrid_df = biogrid_df[mads_interactions]

        # Logging
        logger.info(f'{len(biogrid_df)} MADS-vs-MADS rows retrieved from {len(interactor_list)} UniProt IDs in {self.database}')

        return biogrid_df
    
    def translate_ids(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Translates the BioGRID IDs to UniProt IDs. If several UniProt IDs are
        mapped to a single BioGRID IDs, one is chosen based on being a
        non-deprecated ID and a Swiss-Prot accession. 

        Parameters
        ----------
        df : pd.DataFrame
            BioGRID data frame network with BioGRID ID columns.

        Returns
        -------
        pd.DataFrame
            BioGRID data frame network with UniProt ID coulmns.
        '''
        # Variables used in to_uniprot
        mads_uniprot_ids = [interactor.uniprot_id for interactor in self._iterate_interactors()]
        _, biogrid_to_uniprot = self._get_aliases()
        
        def to_uniprot(biogrid_id: int) -> str:
            '''
            Inner function to apply to each BioGRID ID to map it to a single
            UniProt ID.

            Parameters
            ----------
            biogrid_id : int
                BioGRID ID to map to UniProt IDs

            Returns
            -------
            str
                UniProt ID that best maps the BioGRID ID.
            '''
            
            # Translate BioGRID IDs to UniProt IDs
            uniprot_ids = biogrid_to_uniprot[str(biogrid_id)]

            # Exit if 1-to-1 mapping
            if len(uniprot_ids) == 1:
                return uniprot_ids[0]
            
            # Remove deprecated UniProt IDs
            merged_uniprot_ids = [uniprot_id for uniprot_id in uniprot_ids if uniprot_id in mads_uniprot_ids]

            # Exit if 1-to-1 mapping
            if len(merged_uniprot_ids) == 1:
                return merged_uniprot_ids[0]
            
            # Swiss-Prot or TrEMBL?
            swissprot_uniprot_ids = [uniprot_id for uniprot_id in merged_uniprot_ids if Interactor.unpickle(uniprot_id).uniprot_info['Section'] == 'Swiss-Prot']

            # Exit if a Swiss-Prot accession is found
            if len(swissprot_uniprot_ids) == 1:
                return swissprot_uniprot_ids[0]
            
            # Error to account for BioGRID accessions that map to multiple
            # non-deprecated TrEMBL UniProt IDs
            logger.error(f'BioGRID ID {biogrid_id} cannot be translated properly')
        
        # Apply function to BioGRID ID column
        df['A'] = df['BIOGRID_ID_A'].apply(to_uniprot)
        df['B'] = df['BIOGRID_ID_B'].apply(to_uniprot)   

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
        df['A-B'] = df.apply(lambda x: ' - '.join(sorted([x['A'], x['B']])), axis=1)

        return df


if __name__ == '__main__':
    '''Test class'''
    biogrid = BioGrid()
    #biogrid.map_ids()
    #biogrid.remove_empty_accessions()
    #biogrid.mads_vs_all() 
    df = biogrid.mads_vs_mads()
    df = biogrid.translate_ids(df)
    df = biogrid.standarize(df)
    biogrid.save(df, 'MADS')
    print(df)

    