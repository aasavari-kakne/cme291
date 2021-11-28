from Bio import SeqIO
from Bio.PDB import PDBParser
from collections import OrderedDict
from glob import glob
from itertools import combinations, chain
import numpy as np

from utils import *

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
    def __init__(self, path, parser):
        self.name = get_dir_name(path)
        self.sequence = str(SeqIO.read(path + 'protein.fasta', "fasta").seq)
        self.structure = parser.get_structure("2fat", path + "selected_prediction.pdb")

    def __len__(self):
        return len(self.sequence)

    def get_residues(self):
        model = list(self.structure.get_models())[0]
        chain = list(model.get_chains())[0]
        residues = list(chain.get_residues())
        return residues

    def get_ndarray(self):
        positions = []
        for residue in self.structure:
            for atom in residue.get_atoms():
                p = atom.get_vector().get_array()
                positions.append(p) 
        positions = np.stack(positions, axis=0) if positions else None
        return positions


class Dataset:

    def __init__(self, root):
        self.root = root
        self.parser = PDBParser(PERMISSIVE = True, QUIET = True)
        self.data = self.load_data()
        self.proteins = self.get_proteins()

    def load_data(self):
        data = OrderedDict()
        for gene_path in glob(self.root + '*/'):
            gene_name = get_dir_name(gene_path)
            for protein_path in glob(gene_path + '*/'):
                protein = Protein(protein_path, self.parser)
                data[gene_name] = data.get(gene_name, []) + [protein]
        return data

    def get_proteins(self):
        proteins = []
        for gene in self.data:
            proteins += self.data[gene]
        return proteins