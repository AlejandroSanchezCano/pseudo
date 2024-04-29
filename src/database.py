# Built-in modules
import os
import abc
from collections import defaultdict
from typing import Generator, Literal

# Third-party modules
import pandas as pd
from tqdm import tqdm

# Custom modules
from logger import logger
from interactor import Interactor

# TODO biogrid has no redundancy Uniprot to Biogrid. Does string? NO BUT INTRODUCE WARNING
# TODO ensure that the redundancy from IntAct is correct or that it doesn't bother the mads_vs_all

class Database(abc.ABC):
    
    # Static attributes
    database = None
    uniprot_alias = None
    path = '/home/asanche2/asanche2/ppi/data/'

    @classmethod
    def _get_aliases(cls) -> tuple[dict[str, list[str]], dict[str, list[str]]]: 
        '''
        Open the database file where the database IDs are mapped to external IDs
        (already prefiltered to show UniProt IDs to alleviate some computational
        burden) and return the mapping as a dictionary of list, because I aim to 
        map UniProt IDs to Database IDs and the ID reltionship is many-to-many:
        - 13770 maps to Q9SUP6 and Q0WSJ2 in BioGRID.
        - Q6P996 maps to 3196076 and 116681 in BioGRID.
        Also, because when getting the network, the interactors will be 
        identified with the databases' ID, like in BioGRID or STRING. Therefore,
        we need to map Database IDs with Uniprot ID.
        Returns
        -------
        tuple[dict[str, list[str]], dict[str, list[str]]]
            UniProt ID : list of Database IDs.
            Database ID : list of UniProt IDs
        '''
        # Initialize return variable
        uniprot_to_database = defaultdict(list)
        database_to_uniprot = defaultdict(list)

        # Open database file
        database_file = f'{cls.path}uniprot{cls.uniprot_alias}'
        with open(database_file, 'r') as handle:

            # Parse database file
            for line in tqdm(handle.read().splitlines()):
                database_id, uniprot_id = line.split('\t')

                # Add to alias-mapping dictionary
                uniprot_to_database[uniprot_id] += [database_id]
                database_to_uniprot[database_id] += [uniprot_id]
        
        # Logging
        len_flattened_database = sum(len(mapped) for mapped in uniprot_to_database.values())
        len_flattened_uniprot = sum(len(mapped) for mapped in database_to_uniprot.values())
        logger.info(f'{len(uniprot_to_database)} total UniProt IDs mapped to {len_flattened_database} {cls.database} IDs')
        logger.info(f'{len(database_to_uniprot)} total {cls.database} IDs mapped to {len_flattened_uniprot} UniProt IDs')

        return uniprot_to_database, database_to_uniprot
    
    @classmethod
    def _iterate_interactors(cls) -> Generator[Interactor, None, None]:
        '''
        Iterate over the interactor in a directory.

        Yields
        ------
        Generator[Interactor, None, None]
            Iteractor objects.
        '''
        interactors_path = '/home/asanche2/asanche2/ppi/data/Interactors'
        interactors_pkl = sorted(os.listdir(interactors_path))
        for interactor_pkl in tqdm(interactors_pkl):
            yield Interactor.unpickle(interactor_pkl)

    @classmethod
    def parse_uniprot(cls) -> None:
        '''
        Calls the 'Interactor.parse_uniprot()' method for all Interactors in a 
        directory. It's only meant to be run once.
        '''
        for interactor in cls._iterate_interactors():
            interactor.parse_uniprot()
            interactor.pickle()

    @classmethod
    def map_ids(cls) -> None:
        '''
        Map UniProt IDs to Database ID by accessing a constructed alias 
        dictionary mapping both databases. The Database ID is added to the 
        Interactor object information. 
        Other mapping methods include sequential API requests, with low space
        complexity but bothersome time complexity (2h vs 7min).
        ''' 
        # How many are mapped? (1)
        n_mapped_ids = 0

        # Iterate over UniProt IDs with an alias in Database
        uniprot_to_database, _ = cls._get_aliases()
        for interactor in cls._iterate_interactors():
            if uniprot_to_database.get(interactor.uniprot_id):

                # Update the corresponding method
                match cls.database:
                    case 'BioGRID':
                        interactor.biogrid_id = uniprot_to_database[interactor.uniprot_id]
                    case 'STRING':
                        interactor.string_id = uniprot_to_database[interactor.uniprot_id]
                    case 'IntAct':
                        interactor.intact_id = interactor.uniprot_id
                interactor.pickle()

                # How many are mapped? (2)
                n_mapped_ids += 1
        
        # Logging
        logger.info(f'{n_mapped_ids} MADS UniProt IDs in {cls.database}')
    
    @classmethod
    def in_database(cls) -> Generator[Interactor, None, None]:
        '''
        Uses the 'Interactor().biogrid_id', 'Interactor().string_id', 
        'Interactor().intact_id' or 'Interactor().plappisite_id' previouslly 
        updated by the 'self.map_ids()' method to give the Interactors belonging 
        to a given database.

        Yields
        ------
        Generator[Interactor, None, None]
            Interactors in Database.
        '''
        for interactor in cls._iterate_interactors():
            if cls.database == 'BioGRID' and interactor.biogrid_id:
                yield interactor
            elif cls.database == 'STRING' and interactor.string_id:
                yield interactor
            elif cls.database == 'IntAct' and interactor.intact_id:
                yield interactor
            elif cls.database == 'PlaPPISite' and interactor.plappisite_id:
                yield interactor

    @abc.abstractmethod
    def mads_vs_all(self) -> pd.DataFrame:
        '''
        Accesses database and retrieves the MADS vs. all interactions. 
        Abstract method to be implemented by subclass. 

        Returns
        -------
        pd.DataFrame
            MADS vs. all data frame network.
        '''
        pass

    @abc.abstractmethod
    def mads_vs_mads(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Filters the MADS vs. all data frame to end up with the MADS vs. MADS
        data frame network.
        Abstract method to be implemented by subclass.

        Parameters
        ----------
        df : pd.DataFrame
            MADS vs. all data frame network.

        Returns
        -------
        pd.DataFrame
            MADS vs. MADS data frame network.
        '''
        pass    

    @abc.abstractmethod
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
        pass

    @classmethod
    def save(cls, df: pd.DataFrame, type: Literal['all', 'MADS'] = 'MADS') -> None:
        '''
        Pickles a data frame network to be accessed quickly.

        Parameters
        ----------
        df : pd.DataFrame
            Database data frame network.
        type : Literal['all&', 'MADS'], optional
            Network type can be either MADS vs all ('all') or MADS vs MADS
            ('MADS'), by default 'MADS'.
        '''

        df.to_pickle(f'/home/asanche2/asanche2/ppi/data/Networks/{cls.database}_MADS vs {type}.pkl')