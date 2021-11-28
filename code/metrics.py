from base_model import *
import numpy as np
from procrustes.orthogonal import orthogonal


def overlap_metric(protein1, protein2, data):
    model = Base_model(protein1, protein2)
    data['correct_alignment'] = model.get_correct_alignment()
    data['incorrect_alignment'] = model.get_incorrect_alignment()
    data['overlap_distance'] = model.get_distances()


def procrustes_metric(protein1, protein2, data):
    A, B = protein1.get_ndarray(), protein2.get_ndarray()
    x = orthogonal(A, B)
    d = np.linalg.norm(x.new_a - x.new_b) / max(len(protein1), len(protein2))
    data['procrustes_distance'] = round(d, 2)