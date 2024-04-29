# Built-in modules
import os
import re
from collections import Counter

# Third-party modules
#import deepcoil
import numpy as np
from Bio import PDB
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

# Custom modules
import utils
from logger import logger

class Structure:

    # Genetic code
    three_to_one = {
        'ALA':'A', 'VAL':'V', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F',
        'TYR':'Y', 'TRP':'W', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 
        'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'CYS':'C',
        'GLY':'G', 'PRO':'P'
        }

    def __init__(self, uniprot_id: str):
        
        # Instantiate from 'Structures' folder (from __dict__)
        try:
            structure__dict__ = utils.unpickling(f'/home/asanche2/asanche2/ppi/data/Structures/{uniprot_id.split(".")[0]}.pkl')
            for attribute, value in structure__dict__.items():
                setattr(self, attribute, value)
        
        # Instantiate from 'PDBs' folder (parse PDB)
        except FileNotFoundError:
            self.uniprot_id = uniprot_id.split('.')[0]
            self.path = f'/home/asanche2/asanche2/ppi/data/PDBs/{self.uniprot_id}.pdb'
            structure = PDB.PDBParser().get_structure(id = uniprot_id, file = self.path)
            self.seq3 = '-'.join([residue.resname for residue in structure.get_residues()]) 
            self.seq1 = ''.join([Structure.three_to_one[residue.resname] for residue in structure.get_residues()])
            self.pLDDT = [next(residue.get_atoms()).bfactor for residue in structure.get_residues()]
            #self.secondary_structure = self.dssp()
            #try:
            #    self.domain_architecture_index, self.domain_architecture_residue = self.domain_arquitecture()
            #except IndexError:
            #    logger.info(f'The domain architecture of {self.uniprot_id} could not be determined with DSSP')

    def __repr__(self) -> str:
        return str(self.__dict__)
    
    def range_to_residue(self, range: range) -> str:
        return self.seq1[range[0] : range[-1] + 1]
    
    def kmer_counter(self, k: int = 2):
        return Counter(self.seq1[i : i + k] for i in range(len(self.seq1)+1-k))

    def pLDDT_over_seq(self):
        '''
        Builds a scatterplot from the pLDDT values with respect to the sequence
        index and colors the background according to each pLDDT quality class
        range. At the bottom, draws a rugplot (1-D) of the pLDDT values (as in
        InterPro).
        '''
        # Scatterplot and rugplot
        y = self.pLDDT
        x = list(range(len(y)))
        hue = ['Very high' if i >= 90 else 'Condifent' if i >= 70 else 'Low' if i >= 50 else 'Very low' for i in y]
        plddt_palette = {'Very high': '#0053D6', 'Condifent':'#65CBF3', 'Low':'#FFDB13', 'Very low':'#FF7D45'}
        fig = sns.scatterplot(x=x, y=y, color='black')
        sns.rugplot(x=x, hue = hue, palette = plddt_palette, legend = False, linewidth = 2, expand_margins=True)        

        # Horizontal spans
        ax = plt.gca()
        ax.axhspan(90, 100, alpha = 0.5, color = '#0053D6')
        ax.axhspan(70, 90, alpha = 0.5, color = '#65CBF3')
        ax.axhspan(50, 70, alpha = 0.5, color = '#FFDB13')
        ax.axhspan(0, 50, alpha = 0.5, color = '#FF7D45')
        ax.set_ylim(0, 100)

        # Axis labels
        fig.set_xlabel('Residue index')
        fig.set_ylabel('pLDDT')

        # Show plot
        plt.show() 
    
    def dssp(self) -> str:
        '''
        Infers the secondary structure using DSSP from Bio.

        Returns
        -------
        str
            Secondary structure with DSSP code.
        '''
        dssp = PDB.DSSP(self.structure[0], self.path, dssp='mkdssp')
        _, _, secondary_structure, *_ = zip(*dssp.property_list)
        return''.join(secondary_structure)
        
    def domain_arquitecture(self):

        helixes = list(re.finditer(r'H{4,}', self.secondary_structure))
        domain_architecture_index = {
            'MADS': range(0, helixes[1].start()),
            'MADS helix' : range(*helixes[0].span()), 
            'I': range(*helixes[1].span()),
            'K': range(helixes[2].start(), helixes[3].end()),
            'K 1': range(*helixes[2].span()),
            'K 2': range(*helixes[3].span())
        }

        domain_architecture_residue = {}
        for key, value in domain_architecture_index.items():
            domain_architecture_residue[key] = self.range_to_residue(value)
        
        return domain_architecture_index, domain_architecture_residue

### DEEPCOIL ###

    def _deepcoil(self):
        return utils.unpickling(f'/home/asanche2/asanche2/ppi/data/DeepCoil/{self.uniprot_id}.pkl')
    
    def coiledcoil_propensity(self) -> np.ndarray:
        return self._deepcoil()['cc']
    
    def coiledcoil_predicted(self): # PENSAR QUÉ OUTPUT QUEREMOS 
        propensities = self.coiledcoil_propensity()
        peaks = deepcoil.utils.sharpen_preds(propensities)
        print(peaks)

    def coilcoil_plot(self) -> None:
        results = self._deepcoil()
        deepcoil.utils.plot_preds(results, out_file = 'deepcoil.png')

    def heptad(self, position = 'a'): #TODO: make this pretty once we know what we wanna do with this
        heptad_propensities = self._deepcoil()['hept']
        position_propensities = heptad_propensities[:, 1] if position == 'a' else heptad_propensities[:, 2]
        indeces = np.where(position_propensities > 0.2)[0]
        for i in indeces:
            print(i, self.seq1[i])
        return indeces

    @staticmethod
    def store_structures():
        for pdb in tqdm(os.listdir('/home/asanche2/asanche2/ppi/data/PDBs')):
            structure = Structure(pdb)
            utils.pickling(structure.__dict__, f'/home/asanche2/asanche2/ppi/data/Structures/{pdb.replace("pdb", "pkl")}')
        
    @staticmethod
    def _iterate_structures():
        for structure in tqdm(os.listdir('/home/asanche2/asanche2/ppi/data/Structures')):
            yield Structure(structure)

if __name__ == '__main__':
    Structure.store_structures()


