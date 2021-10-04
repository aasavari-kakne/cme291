# Quantization of difference between 3D protein structures
In this project, we focus on quantitative representation of difference between a pair of protein structures. 


## Report 
1. Add section for alphafold in literature review - explain the need for prediction of protein structures. 

2. What is the biological need for quantization of difference in protein structures ?

3. lDDT measures how well the environment in a reference structure is reproduced in a protein model. It is computed over all pairs of atoms in the reference structure at a distance closer than a predefined threshold R<sub>0</sub>sub> (called inclusion radius), and not belonging to the same residue. These atom pairs define a set of local distances *L*. A distance is considered preserved in the model *M* if it is, within a certain tolerance threshold, the same as the corresponding distance in *L*. If one or both the atoms defining a distance in the set are not present in *M*, the distance is considered non-preserved. For a given threshold, the fraction of preserved distances is calculated. The final lDDT score is the average of four fractions computed using the thresholds 0.5 Å, 1 Å, 2 Å and 4 Å. 


## Presentation
1. Explain the problem statement mathematically. 

2. Explain the need and impact of this project. 

3. Metric of success 

4. Inference 

5. Future goals or tasks

6. 


### Prerequisites

1. Mapping from isoforms(RNA) to proteins(Sequence) is a bijection.  

2. Mapping from gene to isoform is one-to-many. 

3. Any two pair of isoforms transcribed by the same genome will have overlaps among sequences of their proteins. 


### What we need

1. A data directory of gene -> proteins -> .fasta and .pdb files


### Algorithm for technique 1

1. Given sequencing of two proteins (say P<sub>1</sub> and P<sub>2</sub>) which have originated from a single genome, say G.

2. For each protein, say P<sup>(1)</sup>, we have a sequence of amino acids a<sup>(1)</sup><sub>1</sub>, a<sup>(1)</sup><sub>2</sub>, .... a<sup>(1)</sup><sub>n</sub> where n is the length of the protein.

3. We predict structures of the proteins using alphafold2 such that P<sup>(1)</sup><sub>i</sub> are the co-ordinates of a<sup>(1)</sup><sub>i</sub>. 

4. Let us say that proteins P<sup>(1)</sup> and P<sup>(2)</sup> correspond to isoforms I<sup>(1)</sup> and I<sup>(2)</sup>. (Note that : I<sup>(1)</sup> and I<sup>(2)</sup> must be transcribed by the source genome G). 

5. We find the exons shared by I<sup>(1)</sup> and I<sup>(2)</sup>. Then by using the lengths of the exons i.e. number of amino acids in each exon, we find overlapping runs between P<sup>(1)</sup> and P<sup>(2)</sup>. 

6. Thus, we get k overlaps L<sub>1</sub>, L<sub>2</sub>, .... L<sub>k</sub> where each lap L<sub>i</sub> stores starting and ending indices of the overlaps in each protein. I.e. L<sub>i</sub> is represented as tuple of size 4 (i<sub>1</sub>, i<sub>2</sub>, j<sub>1</sub>, j<sub>2</sub>) such that i<sub>2</sub> - i<sub>1</sub> == j<sub>2</sub> - i<sub>1</sub> and a<sup>(1)</sup><sub>i<sub>1</sub> + x </sub> == a<sup>(2)</sup><sub>j<sub>1</sub> + x </sub> for x = 0, 1, ... j<sub>2</sub> - j<sub>1</sub>. 

8. Given an overlapping run for a protein, say a<sup>(1)</sup><sub>i<sub>1</sub></sub> to a<sup>(1)</sup><sub>i<sub>2</sub></sub>, we find distance between each pair of the amino acids to form a symmetric matrix D<sup>(1)</sup><sub>1</sub> such that D<sup>(1)</sup><sub>1</sub>[i][j] = ||a<sup>(1)</sup><sub>i</sub> - a<sup>(1)</sup><sub>j</sub>|| for i, j in {i<sub>1</sub>, .... i<sub>2</sub>}. 

9. Given all such matrices, say D<sup>(1)</sup><sub>1</sub>, D<sup>(1)</sup><sub>2</sub>, ...... D<sup>(1)</sup><sub>m</sub> for P<sub>1</sub> and D<sup>(2)</sup><sub>1</sub>, D<sup>(2)</sup><sub>2</sub>, ...... D<sup>(2)</sup><sub>m</sub> for P<sub>2</sub>, we find sum of norms of differences i.e. ||D<sup>(1)</sup><sub>i</sub> - D<sup>(2)</sup><sub>i</sub>||<sub>F</sub>. 


## Ideas for future 
1. what does it mean to have a different structure ?

2. if two proteins have similar structure, how different structures should they have?

3. can we label proteins to be used in supervised learning -> domains ?

4. If p1 and p2 are similar in structure, will they have same function?

5. mutations in the genome are harmful

6. hypothesis testing

7. significance doesn't imply biological significance

8. size effect -> how large n is 

9. what is in the alphafold error file 

10. sherlock 

11. .gtf files -> code

12. permutation testing / bootstrap

13. paper that compare protein structures or functions

14. procrustes

# 