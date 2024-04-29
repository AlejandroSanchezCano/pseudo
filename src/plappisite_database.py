# Built-in modules
import requests

# Third-party modules
import bs4
import pandas as pd

# Custom modules
from logger import logger
from database import Database

# TODO: STILL PLAPPISITE NEEDS A WAY TO UNEXPLODE THE methods part -> groupby('A-B').agg(list)

class PlaPPISite(Database):

    # Override static attributes
    database = 'PlaPPISite'

    def _soupify(self, uniprot_id: str) -> bs4.BeautifulSoup:
        '''
        PlaPPISite has no API and their downloadable files are not suitable for 
        mapping IDs, so the contents need to be retrieved by parsing the web
        content. Retrieves the web content of a spenific UniProt ID as a
        BeautifulSoup object.

        Parameters
        ----------
        uniprot_id : str
            UniProt ID to search in PlaPPISite.

        Returns
        -------
        bs4.BeautifulSoup
            Soupified web content.
        '''
        url = f'http://zzdlab.com/plappisite/single_idmap.php?protein={uniprot_id}'
        web = requests.get(url).text
        soup = bs4.BeautifulSoup(web, features= 'lxml')
        return soup

    def _get_table(self, soup: bs4.BeautifulSoup) -> pd.DataFrame:
        '''
        Parses the web content of an accession and retrieves the PPI table.

        Parameters
        ----------
        soup : bs4.BeautifulSoup
            Accession's web content.

        Returns
        -------
        pd.DataFrame
            Interaction table.
        '''
        table = soup.find('div', attrs = {'id':'container_table'})
        columns = [th.text for th in table.find_all('th')]
        tds = [td.text for td in table.find_all('td')]
        rows = [tds[i : i + len(columns)] for i in range(0, len(tds), len(columns))]

        return pd.DataFrame(rows, columns = columns)

    def map_ids(self) -> None:
        '''
        Map UniProt IDs to PlaPPISite IDs to know the UniProt IDs for which 
        PlaPPISite has interaction data. PlaPPISite uses UniProt IDs as main IDs.
        '''
        # How many are mapped? (1)
        n_mapped_ids = 0

        # Retrieve PPI table of all MADS proteins   
        for interactor in self._iterate_interactors():
            soup = self._soupify(interactor.uniprot_id)
            table = self._get_table(soup)
            
            # Mapped ID when PPI table is not empty
            if not table.empty:
                
                # Update Interactor method
                interactor.plappisite_id = interactor.uniprot_id
                interactor.pickle()

                # How many are mapped? (2)
                n_mapped_ids += 1

        # Logging
        logger.info(f'{n_mapped_ids} MADS UniProt IDs in {self.database}')

    def remove_predicted(self) -> None:
        '''
        PlaPPISite contains both experimental and predicted interactions. Only
        the experimental PPIs are of interest because STRING is arguably the
        best source of predicted PPIs, so PlaPPISite predicted PPIs will not
        be likely used.
        '''
        # How many are removed? (1)
        n_removed = 0
        
        # Retrieved PPI table of already mapped proteins
        for interactor in self.in_database():
            soup = self._soupify(interactor.uniprot_id)
            table = self._get_table(soup)
            non_predicted_table = table[table['PPI source'].apply(lambda x: x not in ['Predicted', 'prediction'])]
            
            # Mapped ID when accession has experimental PPIs
            if non_predicted_table.empty:

                # Update Interactor method
                interactor.plappisite_id = ''
                interactor.pickle()

                # How many are removed? (2)
                n_removed += 1
                
        # Logging
        logger.info(f'{n_removed} MADS UniProt IDs were remove because all its interactions were predicted interactions')

    def mads_vs_all(self) -> pd.DataFrame:
        '''
        Accesses web content for all the mapped UniProt IDs and retrieves the 
        MADS vs. all interactions. 

        Returns
        -------
        pd.DataFrame
            MADS vs. all data frame network.
        '''
        # Initialize MADS vs. all data frame
        df = pd.DataFrame()

        # Loop over mapped IDs.
        interactors = list(self.in_database())
        for interactor in interactors:
            
            # Find experimental PPIs
            soup = self._soupify(interactor.uniprot_id)
            table = self._get_table(soup)
            non_predicted_table = table[table['PPI source'].apply(lambda x: x not in ['Predicted', 'prediction'])]
            
            # Concatenate with main data frame
            df = pd.concat([df, non_predicted_table])
        
        # Logging
        logger.info(f'{len(df)} MADS-vs-all rows retrieved from {len(interactors)} UniProt IDs in {self.database}')

        return df
    
    def mads_vs_mads(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Filters the MADS vs. all data frame to end up with the MADS vs. MADS
        data frame network.

        Parameters
        ----------
        df : pd.DataFrame
            MADS vs. all data frame network.

        Returns
        -------
        pd.DataFrame
            MADS vs. MADS data frame network.
        '''
        # MADS PlaPPISite IDs
        interactor_list = [interactor for interactor in self.in_database()]
        mads_plappisite_ids = set([interactor.plappisite_id for interactor in interactor_list])

        # Filter MADS vs. all network
        df[['A', 'B']] = df['PPI'].str.split(' - ', expand = True)
        df = df.loc[df['A'].apply(lambda x: x in mads_plappisite_ids)]
        df = df.loc[df['B'].apply(lambda x: x in mads_plappisite_ids)]

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
        df = df.rename(columns = {'PPI' : 'A-B'})

        return df
            
if __name__ == '__main__':
    plappisite = PlaPPISite()
    #plappisite.map_ids()
    #plappisite.remove_predicted()
    df = plappisite.mads_vs_all()
    df = plappisite.mads_vs_mads(df)
    df = plappisite.standarize(df)
    plappisite.save(df, 'MADS')
    print(df)