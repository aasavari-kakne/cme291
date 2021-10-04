import numpy as np
from protein import *
from alignment import *

class Base_model:
    

    def __init__(self, protein1_path, protein2_path):
        # read proteins
        self.protein1 = Protein(protein1_path) 
        self.protein2 = Protein(protein2_path)
        # read their structures
        self.res1 = self.protein1.get_residues()
        self.res2 = self.protein2.get_residues()
        # perform alignment
        self.alignment = Alignment(self.protein1.sequence, self.protein2.sequence)
        # get overlaps from alignment
        self.overlaps = self.alignment.get_overlaps()
        # get mapping from indices of alignment to indices of original amino acid sequence
        self.map1, self.map2 = self.alignment.get_index_maps()
        

    def get_positions_for_overlap(self, residues, start, end):
        positions = []
        for i in range(start, end + 1):
            for atom in residues[i].get_atoms():
                p = atom.get_vector().get_array()
                positions.append(p) 
        positions = np.stack(positions, axis=0) if positions else None
        return positions
    

    def get_dist_mat(self, p):
        z = p[:, np.newaxis, :] - p[np.newaxis, :]
        D = np.linalg.norm(z, axis=2)
        return D

    def round(self, num):
        return float(np.round(num, 2))


    def get_alignment_percentage(self):
        percent_alignment = 100.0 * self.alignment.get_score() / max(len(self.protein1), len(self.protein2))
        return self.round(percent_alignment)
    

    def get_distances(self):
        distance = 0
        matrices = []
        for s, e in self.overlaps:
            p1 = self.get_positions_for_overlap(self.res1, self.map1[s], self.map1[e]) # n x 3
            p2 = self.get_positions_for_overlap(self.res2, self.map2[s], self.map2[e]) # n x 3
            D1 = self.get_dist_mat(p1)
            D2 = self.get_dist_mat(p2)
            distance += np.linalg.norm(D1 - D2)
            matrices += (D1, D2),
        return matrices, self.round(distance)
