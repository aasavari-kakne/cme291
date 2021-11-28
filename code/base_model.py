import numpy as np
from alignment import *

class Base_model:
    """

        implements baseline model for computing distance between two proteins


        attributes :
        ------------

            protein1 : Protein to store sequence and structure of first protein

            protein2 : Protein to store sequence and structure of second protein

            res1 : list of residues from first protein

            res2 : list of residues from second protein

            alignment : Alignment to store and process alignment between two 
                        proteins

            overlaps : list of tuples which store starting and ending indices
                       of valid overlaps between two proteins 

            map1 : list to map indices of alignment to original sequence of 
                   first protein

            map2 : list to map indices of alignment to original sequence of 
                   second protein


        methods :
        ---------

            get_correct_alignment : return number of correctly aligned residue 
                                    pairs

            get_incorrect_alignment : return number of incorrectly aligned residue 
                                    pairs

            get_positions_for_overlap : returns ndarray containing 3D positions 
                                        of all atoms in a given overlap


            round : rounds floats upto 2 decimal places


            get_dist_mat : returns 2D matrix D where D[i, j] = || p[i] - p[j] ||
                            and p[i] stores 3D position of atom i.
            
            get_distances : returns a list of tuple of distance matrices (D1, D2)
                            such that each tuple corresponds to an overlap 
                            between the amino acid sequence of the two proteins.
                            Also, returns distance values for defined metrics.

    """
    def __init__(self, protein1, protein2):
        # read proteins
        self.protein1 = protein1
        self.protein2 = protein2
        # read their structures
        self.res1 = self.protein1.get_residues()
        self.res2 = self.protein2.get_residues()
        # perform alignment
        self.alignment = Alignment(self.protein1.sequence, self.protein2.sequence)
        # get overlaps from alignment
        self.overlaps = self.alignment.get_overlaps()
        # get mapping from indices of alignment to indices of original sequence
        self.map1, self.map2 = self.alignment.get_index_maps()

    def get_correct_alignment(self):
        return self.alignment.get_score()

    def get_incorrect_alignment(self):
        return max(len(self.protein1), len(self.protein2)) - self.alignment.get_score()

    def get_positions_for_overlap(self, residues, start, end):
        positions = []
        for i in range(start, end + 1):
            for atom in residues[i].get_atoms():
                p = atom.get_vector().get_array()
                positions.append(p) 
        positions = np.stack(positions, axis=0) if positions else None
        return positions

    def round(self, num):
        return float(np.round(num, 2))
    
    def get_dist_mat(self, p):
        z = p[:, np.newaxis, :] - p[np.newaxis, :]
        D = np.linalg.norm(z, axis=2)
        return D

    def get_matrices(self):
        matrices = []
        for s, e in self.overlaps:
            p1 = self.get_positions_for_overlap(self.res1, self.map1[s], self.map1[e])
            p2 = self.get_positions_for_overlap(self.res2, self.map2[s], self.map2[e])
            if p1.shape == p2.shape:
                matrices += (p1, p2),
        return matrices

    def get_distances(self):
        distance = 0
        for p1, p2 in self.get_matrices():
                D1 = self.get_dist_mat(p1)
                D2 = self.get_dist_mat(p2)
                d = np.linalg.norm(D1 - D2) / p1.shape[0]
                distance += d
        return self.round(distance)

