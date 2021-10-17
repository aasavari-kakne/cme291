from Bio import pairwise2

class Alignment:
	"""

		class to generate and store alignment between two protein sequences
		

		attributes :

			alignment : Bio.pairwise2.Alignment to store the alignment between 
						two proteins

		methods :

			__len__ : returns the length of alignment 

			get_score : returns the number of correctly aligned amino acid pairs 
						between two proteins

			get_index_map : creates a mapping from indices of aligned sequence 
							to the original sequence 
							>>> get_index_map('AC-T') 
							[0, 1, -1, 2]

			get_index_maps : returns index maps for the aligned sequences 
							 which are stored in the alignment attribute

			get_alignment_mask : returns a mask with boolean values, where
								 each value indicates if there is a correct 
								 alignment at this index or not

			get_overlaps : returns a list of tuples where each tuple stores 
						   two indices. First index indicates starting of a valid
						   overlap between two proteins. Second index indicates 
						   ending of the overlap (both indices inclusive)


	"""
	def __init__(self, seq1, seq2):
		self.alignment = max(pairwise2.align.globalxx(seq1, seq2))


	def __len__(self):
		return len(self.alignment.seqA)


	def get_score(self):
		return self.alignment.score
		

	def get_index_map(self, seq):
		n = len(seq)
		start, mapping = 0, [-1]*n
		for i in range(n):
			if seq[i].isalpha():
				mapping[i] = start
				start += 1
		return mapping
	

	def get_index_maps(self):
		return [self.get_index_map(seq) for seq in [self.alignment.seqA, self.alignment.seqB]]


	def get_alignment_mask(self):
		n = len(self)
		mask = [self.alignment.seqA[i] == self.alignment.seqB[i] for i in range(n)]
		return mask


	def get_overlaps(self):
		n, start, overlaps = len(self), 0, []
		mask = self.get_alignment_mask()
		while start < n:
			# find starting index of an overlap
			while start < n and mask[start] == False:
				start += 1
			if start == n:
				break
			# find index after the ending of the overlap
			end = start
			while end < n and mask[end] == True:
				end += 1
			# update list of overlaps 
			overlaps += (start, end - 1),
			# update starting point
			start = end
		return overlaps

