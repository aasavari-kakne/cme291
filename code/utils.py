import numpy as np

def get_positions(data):
	"""
		given an instance of Bio.PDB.Structure.Structure, 
		collects and returns positions of all the atoms in the structure

		args:
			data : instance of class Bio.PDB.Structure.Structure 

		returns:
			positions : ndarray of size (n, 3),
						where n is number of atoms in the structure

	"""
	positions = []
	for model in data.get_models():
		for chain in model.get_chains():
			for residue in chain.get_residues():
				for atom in residue.get_atoms():
					p = atom.get_vector().get_array()
					positions.append(p) 
	positions = np.stack(positions, axis=0) if positions else None
	return positions


def get_dist_mat(p):
	"""
		given a position array, returns the corresponding distance matrix 
		
		args:
			p : ndarray of size (n, 3) where n is number of atoms 
			
		returns:
			D : ndarray of size (n, n),
				where D[i, j] = || p[i] - p[j] || for i, j = 0, 1, ... n-1
			
	"""
	z = p[:, np.newaxis, :] - p[np.newaxis, :]
	D = np.linalg.norm(z, axis=2)
	return D


def get_index_map(seq):
	"""
		maps indices of the aligned sequence to it's corresponding index in the 
		original sequence

		args:
			seq : aligned sequence of size n 
				  where each value is either an amino acid or '-' 

		returns:
			mapping : list of size n where each item either holds 
			a valid index in the original sequence 
			or
			-1 which corresponds to '-' in seq
	"""
	n = len(seq)
	start, mapping = 0, [-1]*n
	for i in range(n):
		if seq[i].isalpha():
			mapping[i] = start
			start += 1
	return mapping


def get_overlaps(alignment):
	"""
		returns list of tuples of size 2 e.g. (start, end)
		where residues [s:e+1] indicate valid alignment 
	"""
	n = len(alignment.seqA)
	mask = [int(alignment.seqA[i] == alignment.seqB[i]) for i in range(n)]
	overlaps = []
	start = 0
	while start < n:
		# find starting index of an overlap
		while start < n and mask[start] == 0:
			start += 1
		if start == n:
			break
		# find index after the ending of the overlap
		end = start
		while end < n and mask[end]:
			end += 1
		overlaps += (start, end - 1),
		start = end
	return overlaps

def get_residues(data):
    model = list(data.get_models())[0]
    chain = list(model.get_chains())[0]
    residues = list(chain.get_residues())
    return residues
    

