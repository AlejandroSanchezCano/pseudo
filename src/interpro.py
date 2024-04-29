# Built-in modules
import os
import math
import time
import subprocess
from typing import Any

# Third-party modules
import requests
from tqdm import tqdm

# Custom modules
import utils
from logger import logger
from interactor import Interactor

class InterPro:

    # Static attributes
    url = "https://www.ebi.ac.uk/interpro/api"
    interpro2domain = {
        'IPR036879' : 'MADS-box',
        'IPR002487' : 'K-box'
    }

    def __init__(self, accession: str):
        # Modified upon instantiation
        self.accession = accession
        self.source_database = self.__get_source_database()

        # Modified with 'get_metadata()' method
        self.name = ''
        self.number_of = {
            'proteins': 0,
            'alphafolds': 0
        }

        # Modified with 'get_uniprot()' method
        self.interactors = []

    def __get_source_database(self) -> str:
        '''
        Parses the Interpro accession ID to find what database it belongs to
        within InterPro (CDD, Profile, Panther, PFAM ...). It returns the 
        corresponding string used in API requests.

        Returns
        -------
        str
            InterPro
        '''
        if self.accession.startswith('IPR'):
            return 'interpro'
        elif self.accession.startswith('cd'):
            return 'cdd'
        elif self.accession.startswith('G3DSA'):
            return 'cathgene3d'
        elif self.accession.startswith('P'):
            return 'profile'
    
    def __repr__(self) -> str:
        '''
        Prints instance attributes for debugging purposes.

        Returns
        -------
        str
            Instance attribute mapping.
        '''
        return str(self.__dict__)

    def __request(self, url: str) -> dict[str, Any]:
        '''
        Given an InterPro API url, performs a request and returns the response
        as a python-interactable JSON. 

        Parameters
        ----------
        url : str
            InterPro API url.

        Returns
        -------
        dict[str, Any]
            Parsed JSON response.
        '''
        request = requests.get(url)
        return request.json()
    
    def get_metadata(self) -> None:
        '''
        Performs an InterPro API request to retrieve the name of the accession
        and the number of UniProt and AlphaFold accessions it contains. Updates
        'self.name' and 'self.number_of' attributes.
        '''
        # Request to InterPro API
        url = f'{InterPro.url}/entry/interpro/{self.accession}'
        json = self.__request(url)
        
        # Access metadata in API response
        self.name = json['metadata']['name']['name']
        self.number_of['proteins'] = json['metadata']['counters']['proteins']
        self.number_of['alphafolds'] = json['metadata']['counters']['structural_models']['alphafold']

        # Logging
        logger.info(f'{self.name=}')
        logger.info(f'{self.number_of=}')
    
    ''' # We might not even need this cause it's already in metadata
    def __get_total_batches(self, batch_size: int = 200) -> int:
        #InterPro API manages large API requests by pagination and providing the
        #url to the next page in the previous page, 

        # Separate because otherwise we cannot have the information of the batches in tqdm
        url = f'{InterPro.url}/protein/uniprot/entry/{self.source_database}/{self.accession}?page_size=1'
        json = self.__request(url)
        n_uniprot_entries = json['count']
        return math.ceil(n_uniprot_entries/batch_size)
    '''

    def __get_domains(self, json_result: dict[str, Any]) -> list[tuple[int, int]]: # Interesting -> A0A061FZ27 has two K-box domains!!
        '''
        Inspect the JSON info of the InterPro protein to retrieve the domain
        bounds. 

        Parameters
        ----------
        json_result : dict[str, Any]
            Protein JSON.

        Returns
        -------
        list[tuple[int, int]]
            List of domain bounds as tuples of (start, end).
        '''
        # Initialize return variable
        domains = []

        # Assert that there is only 1 'entry' and 'fragment' per protein
        uniprot_id = json_result['metadata']['accession']
        assert len(json_result['entries']) == 1, f'More than one "entries" in {uniprot_id}'
        assert len(json_result['entries'][0]['entry_protein_locations'][0]['fragments']) == 1, f'More than one "fragments" in {uniprot_id}'

        # Find start and end of domains
        for entry_protein_location in json_result['entries'][0]['entry_protein_locations']:
            start = entry_protein_location['fragments'][0]['start'] - 1  # InterPro is 1-based
            end = entry_protein_location['fragments'][0]['end'] - 1
            domains += [(start, end)]

        return domains

    def get_uniprot(self, batch_size: int = 200) -> None:
        '''
        Use InterPro API to retrieve necessary the UniProt IDs, taxons and 
        domains of the proteins belonging to the self.accession InterPro ID.

        Parameters
        ----------
        batch_size : int, optional
            Batch size, by default 200
        '''

        # InterPro API URL
        url = f'{InterPro.url}/protein/uniprot/entry/{self.source_database}/{self.accession}?page_size={batch_size}'

        # Manage API pagination
        total_batches = math.ceil(self.number_of['proteins']/batch_size)
        with tqdm(total = total_batches) as pbar:
            while url:
                
                # Page response JSON
                json = self.__request(url)
                for result in json['results']:
                    
                    # Access protein info
                    uniprot_id = result['metadata']['accession']
                    taxon_id = result['metadata']['source_organism']['taxId']
                    domains = self.__get_domains(result)

                    # Initialize Interactor object with InterPro info
                    interactor = Interactor(
                        uniprot_id = uniprot_id, 
                        taxon_id = taxon_id,
                        domains = domains
                        )
                    
                    self.interactors.append(interactor)

                # Prepare for next batch
                url = json['next']
                pbar.update(1)
                time.sleep(0.2)

    def __and__(self, other) -> list[Interactor]:
        '''
        Manages intersection behaviour between the Interactor of two InterPro
        instances. It is used to find the MIKC proteins from the MADS-containing
        proteins and the K-box-containing proteins.
        If slow, try "awk 'NR==FNR{arr[$0];next} $0 in arr' %s.txt %s.txt | wc -l"

        Returns
        -------
        list[Interactor]
            List of Interactor objects corresponding to the common ones shared
            by two InterPro instances.
        '''
        return list(set(self.interactors) & set(other.interactors))

    @staticmethod
    def save(interactors: list[Interactor]) -> None:
        '''
        Pickle the __dict__ of a list of Interactor objects (i.e. MIKC).

        Parameters
        ----------
        interactors : list[Interactor]
            List of Interactor objects to pickle.
        '''
        for interactor in interactors:
            utils.pickling(data = interactor.__dict__, path = f'/home/asanche2/asanche2/ppi/data/Interactors/{interactor.uniprot_id}.pkl')

    def get_structures(self, interactors: list[Interactor]) -> None: #14106 out of 16990
        '''
        Download AlphaFold structures from the AlphaFold database of the list of
        MIKC proteins.

        Parameters
        ----------
        interactors : list[Interactor]
            List of Interactor objects to download structures of.
        '''
        for interactor in tqdm(interactors):
            # Download response from AlphaFold database
            structure_path = f'/home/asanche2/asanche2/ppi/data/Structures/{interactor.uniprot_id}.pdb'
            cmd = f'curl "https://alphafold.ebi.ac.uk/files/AF-{interactor.uniprot_id}-F1-model_v4.pdb"'
            response = subprocess.run(cmd, capture_output = True, text = True, shell = True).stdout

            # Only save the non-empty responses (> 200 characters)
            if len(response) > 200:
                with open(structure_path, 'w') as handle:
                    handle.write(response)

        # Logging
        n_pdb_files = len(os.listdir(f'/home/asanche2/asanche2/ppi/data/Structures/'))
        logger.info(f'{n_pdb_files} AlphaFold files have been downloaded out of {len(interactors)}')

if __name__ == '__main__':
    '''Test class'''
    # MADS proteins
    kbox = InterPro('IPR002487')
    kbox.get_metadata()
    kbox.get_uniprot()

    # K-box proteins
    mads = InterPro('IPR036879')
    mads.get_metadata()
    mads.get_uniprot()
    
    # MIKC proteins
    mikc = kbox & mads
    InterPro.save(mikc)
    InterPro.get_structures(mikc)