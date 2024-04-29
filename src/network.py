# Built-in modules
from collections import defaultdict

# Third-party modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Custom modules
import utils
from logger import logger
from structure import Structure
from interactor import Interactor
from interaction import Interaction

class Network:

    def __init__(self, database = 'STRING', type = 'MADS vs MADS', species = 'all'):
        self.database = database
        self.path = f'/home/asanche2/asanche2/ppi/data/Networks/{database}_{type}.pkl'
        self.df = utils.unpickling(self.path) # EL HOMOGENIZE SE HACE EN CADA DB 
        self.df = self.df.drop_duplicates('A-B') # CAREFUL CAUSE WE LOSE INFO
        self.df = self.add_species()
        #self.df['A-B'] = self.df['A-B'].apply(Interaction)
        self.df = self.species(species) if species != 'all' else self.df

    def __len__(self) -> int:
        return len(self.df)
    
    def interactors(self) -> set:
        return list(set(self.df.A) | set(self.df.B))

    def interaction_density(self) -> int:
        ppi = len(self)
        all_interactors = len(self.interactors())
        posible_interactions = sum(range(all_interactors))
        return ppi/posible_interactions
    
    def add_species(self):
        self.df['Species A'] = self.df['A'].apply(lambda x : Interactor.unpickle(x).taxon_id)
        self.df['Species B'] = self.df['B'].apply(lambda x : Interactor.unpickle(x).taxon_id)

        return self.df

    def species(self, *species: int):
        species_A = self.df['Species A'].isin(species)
        species_B = self.df['Species B'].isin(species)
        return self.df[species_A & species_B]
    
    def plot_per_taxa(self):
        taxa_freq = defaultdict(int)

        for interactor in self.interactors():
            interactor = Interactor.unpickle(interactor)
            taxa_freq[str(interactor.taxon_id)] += 1  # taxon_id_to_name when multitax works
        
        taxa_freq = dict(sorted(taxa_freq.items(), key = lambda x: x[1], reverse=True))
        
        plt.bar(*zip(*taxa_freq.items()))
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.xlabel('Species')
        plt.ylabel('Number of UniProt accessions')
        plt.title('')
        plt.show()
    
    def positive_pairs(self, order: bool = False) -> set[str]:
        # We don't care about order 
        if not order:
            return set(self.df['A-B'])
        
        # We care about order
        positive_pairs = []
        for pair in self.df['A-B']:
            A, B = pair.split('-')
            AB = A + '-' + B
            BA = B + '-' + A
            
            positive_pairs += [AB, BA]
        
        return set(positive_pairs)

    def negative_pairs(self, order: bool = False) -> set[str]:
        negative_pairs = []
        positive_pairs = self.positive_pairs(order = order)
        interactors = self.interactors()
        for pointer_A in range(len(interactors)):
            for pointer_B in range(pointer_A, len(interactors)):
                A = interactors[pointer_A]
                B = interactors[pointer_B]
                AB = A + '-' + B
                BA = B + '-' + A

                if not order:
                    negative_pairs += [] if AB in positive_pairs or BA in positive_pairs else [AB]
                else:
                    negative_pairs += [] if AB in positive_pairs else [AB, BA]
        
        return set(negative_pairs)

if __name__ == '__main__':
    n = Network(database = 'IntAct', type = 'MADS vs MADS', species = 'all')
    print(len(n.interactors()))
