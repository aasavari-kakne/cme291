from Bio import pairwise2

class Alignment:
    
    def __init__(self, seq1, seq2):
        self.alignment = max(pairwise2.align.globalxx(seq1, seq2))
        
    def get_overlaps(self):
        n = len(self.alignment.seqA)
        mask = [int(self.alignment.seqA[i] == self.alignment.seqB[i]) for i in range(n)]
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
            while end < n and mask[end] == 1:
                end += 1
            overlaps += (start, end - 1),
            start = end
        return overlaps
    
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

    def get_score(self):
        return self.alignment.score