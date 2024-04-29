import deepcoil
import utils
import os
from tqdm import tqdm
#from deepcoil import DeepCoil
#from deepcoil.utils import plot_preds
#from deepcoil.utils import sharpen_preds 
from Bio import SeqIO
import numpy as np
from structure import Structure
#dc = DeepCoil(use_gpu=True)

class CoiledCoilPrediction():

    def __init__(self):
        self.deepcoil = deepcoil.DeepCoil(use_gpu = False)

    def _predict(self, structure):
        fasta = {structure.uniprot_id : structure.seq1}
        results = self.deepcoil.predict(fasta)[structure.uniprot_id]
        
        return results

    def _pickle(self, results, uniprot_id) -> None:
        utils.pickling(results, f'/home/asanche2/asanche2/ppi/data/DeepCoil/{uniprot_id}.pkl')

    def _iterate_structures(self):
        structures_path = '/home/asanche2/asanche2/ppi/data/Structures'
        structures_pkl = sorted(os.listdir(structures_path))
        for structure_pkl in tqdm(structures_pkl):
            yield Structure(structure_pkl)
    
    def predict_all(self):
        for structure in self._iterate_structures():
            results = self._predict(structure)
            self._pickle(results, structure.uniprot_id)

if __name__ == '__main__':
    CoiledCoilPrediction().predict_all()