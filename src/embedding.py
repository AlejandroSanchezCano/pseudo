# Built-in modules
from collections import Counter

# Third-party modules
import numpy as np
import pandas as pd
from sklearn.preprocessing import MaxAbsScaler

# Custom modules
from logger import logger
from structure import Structure
from network import Network
from interaction import Interaction

class Embedding:

    def __init__(self, network):
        self.network = network
        self.df = network.df
        self.interactors = network.interactors()

    def uniprotid_to_kmer(self):
        dictionary = {}
        for uniprot_id in self.interactors:
            dictionary[uniprot_id] = Structure(uniprot_id).kmer_counter()

        return dictionary
    
    def maxabs_scaler(self, embedding):
        features = embedding.columns
        proteins = embedding.index
        scaled_embedding = MaxAbsScaler().fit_transform(embedding) 

        return pd.DataFrame(scaled_embedding, index = proteins, columns = features)
    
    def length_scaler(self, embedding):
        pair_length = embedding.sum(axis = 1)  
        scaled_embedding = embedding.div(pair_length, axis = 0)

        return scaled_embedding

    def kmer_sum(self):
        interaction_to_kmer = {}
        uniprotid_to_kmer = self.uniprotid_to_kmer()
        total_pairs = self.network.positive_pairs(0) | self.network.negative_pairs(0)

        for interaction in total_pairs:
            interaction = Interaction(interaction)
            kmer_a = uniprotid_to_kmer[interaction.A.uniprot_id]
            kmer_b = uniprotid_to_kmer[interaction.B.uniprot_id]
            kmers_a_b = sum([kmer_a, kmer_b], Counter())
            interaction_to_kmer[interaction] = kmers_a_b
        
        return pd.DataFrame(interaction_to_kmer).T.replace(np.nan, 0)
    
    def kmer_concat(self):
        interaction_to_kmer = {}
        uniprotid_to_kmer = self.uniprotid_to_kmer()
        total_pairs = self.network.positive_pairs(1) | self.network.negative_pairs(1)
        
        
        for interaction in total_pairs:
            interaction = Interaction(interaction)
            kmer_a = uniprotid_to_kmer[interaction.A.uniprot_id]
            kmer_b = uniprotid_to_kmer[interaction.B.uniprot_id]

            # A-B
            kmer_a_A = {f'{kmer}_A':kmer_a[kmer] for kmer in kmer_a}
            kmer_b_B = {f'{kmer}_B':kmer_b[kmer] for kmer in kmer_b}
            kmer_concat = kmer_a_A | kmer_b_B
            interaction = Interaction(interaction.A.uniprot_id + '-' + interaction.B.uniprot_id)
            interaction_to_kmer[interaction] = kmer_concat

            # B-A
            kmer_a_B = {f'{kmer}_B':kmer_a[kmer] for kmer in kmer_a}
            kmer_b_A = {f'{kmer}_A':kmer_b[kmer] for kmer in kmer_b}
            kmer_concat = kmer_a_B | kmer_b_A
            interaction = Interaction(interaction.B.uniprot_id + '-' + interaction.A.uniprot_id)
            interaction_to_kmer[interaction] = kmer_concat

        return pd.DataFrame(interaction_to_kmer).T.replace(np.nan, 0)
    
    def multiindex(self, df: pd.DataFrame) -> pd.DataFrame:
        # Initialize variables
        pair_to_index = {} 
        counter_index = 0

        # Map the A-B and B-A pairs to an index number 
        for interaction in df.index:
            if interaction.B_A not in pair_to_index:
                pair_to_index[interaction.A_B] = counter_index
                counter_index += 1
            else:
                pair_to_index[interaction.A_B] = pair_to_index[interaction.B_A]

        # Interaction object mapped to index number
        interaction_to_index = {Interaction(k):v for k,v in pair_to_index.items()}

        # Multiindex
        multiindex = pd.DataFrame(interaction_to_index.items(), columns = ['Interaction', 'Index'])
        multiindex = multiindex.sort_values('Index').loc[:, ['Index', 'Interaction']]
        multiindex = pd.MultiIndex.from_frame(multiindex)

        # Add to embedding
        df_multiindex = pd.DataFrame(df.to_numpy(), index = multiindex, columns = df.columns)

        return df_multiindex
    
    def response(self, embedding):
        positive_pairs = self.network.positive_pairs(order = True)
        interactions = embedding.index.get_level_values('Interaction').to_series()
        interactions = interactions.apply(lambda x: x.A_B)
        response = interactions.isin(positive_pairs).map(int)
        response.index = embedding.index
        return response

if __name__ == '__main__':
    intact_network = Network(database = 'IntAct', type = 'MADS vs MADS', species = 'all')

    e = Embedding(intact_network)
    kmer_concat = e.kmer_concat()
    print(kmer_concat)
    df = e.multiindex(kmer_concat)
    e.response(df)

    print(df)
    #kmer_sum = e.length_scaler(kmer_sum)
    #kmer_sum = e.maxabs_scaler(kmer_sum) * 100
    #y = e.response(kmer_sum)

    #kmer_concat = e.kmer_concat()
    #kmer_concat = e.length_scaler(kmer_concat)
    #kmer_concat = e.maxabs_scaler(kmer_concat) * 100

    
    #from machine_learning import MachineLearning
    #X = kmer_sum
    #y = e.response(kmer_sum)
    #print(MachineLearning(X, y).lazy())