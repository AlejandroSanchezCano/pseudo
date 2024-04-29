from interactor import Interactor

class Interaction:

    def __init__(self, A_B: str):
        self.A_B = A_B
        self.A, self.B = [Interactor.unpickle(interactor) for interactor in A_B.split('-')]
        self.B_A = self.B.uniprot_id + '-' + self.A.uniprot_id
        self.species = self.A.taxon_id if self.A.taxon_id == self.B.taxon_id else 'interspecies'

    def __repr__(self):
        return self.A_B
    
    def __hash__(self):
        return hash(self.A_B)
    
    def __eq__(self, other):
        return self.A_B == other.A_B
        
if __name__ == '__main__':
    i = Interaction('F4KCU5-P29382')
    j = Interaction('P29382-F4KCU5')
    print(i == j)
    print(i.B_A)