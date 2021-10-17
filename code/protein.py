from Bio import SeqIO
from Bio.PDB import PDBParser

class Protein:
    """
        reads and stores attributes of a protein 

        attributes :

            sequence : str to store the sequence of amino acids in the protein

            structure : Bio.PDB.Structure.Structure to store the 3D structure 
            of the protein, as predicted by the alphafold model

        methods :

            __len__ : returns number of amino acids in the protein

            get_residues : returns list of residues for the first chain 
            of the first model of the protein structure
    """
    def __init__(self, path):
        # read sequence of amino acids
        self.sequence = str(SeqIO.read(path + 'protein.fasta', "fasta").seq)
        # read 3D structure of protein
        parser = PDBParser(PERMISSIVE = True, QUIET = True)
        self.structure = parser.get_structure("2fat", path + "selected_prediction.pdb")


    def __len__(self):
        return len(self.sequence)


    def get_residues(self):
        model = list(self.structure.get_models())[0]
        chain = list(model.get_chains())[0]
        residues = list(chain.get_residues())
        return residues
