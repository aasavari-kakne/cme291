# Quantization of difference between 3D protein structures
In this project, we focus on quantitative representation of difference between a pair of protein structures. Please follow the high-level stucture and diagrams of code as well as metric formulations with explanation from the presentation here - https://docs.google.com/presentation/d/1niliTaGCRhdtE83dPnewiufjewYER91gT8Rn_qqWNoQ/edit?usp=sharing. 


## Literature review:
1. Please find and closely follow the Alphafold and Alphafold companion paper closely in papers folder. 
2. **Important note regarding the evaluation metric of Alphafold**: lDDT measures how well the environment in a reference structure is reproduced in a protein model. It is computed over all pairs of atoms in the reference structure at a distance closer than a predefined threshold R<sub>0</sub>sub> (called inclusion radius), and not belonging to the same residue. These atom pairs define a set of local distances *L*. A distance is considered preserved in the model *M* if it is, within a certain tolerance threshold, the same as the corresponding distance in *L*. If one or both the atoms defining a distance in the set are not present in *M*, the distance is considered non-preserved. For a given threshold, the fraction of preserved distances is calculated. The final lDDT score is the average of four fractions computed using the thresholds 0.5 Å, 1 Å, 2 Å and 4 Å. 


### Important assumptions of the project
1. Mapping from isoforms(RNA) to proteins(Sequence) is a bijection.  
2. Mapping from gene to isoform is one-to-many. 
3. Any two pair of isoforms transcribed by the same genome will have overlaps among sequences of their proteins. 


### Input data 

1. We need a data directory of root -> gene -> proteins -> .fasta and .pdb files


## Minutes of meetings (High level Ideas for future)
1. what does it mean to have a different structure -> map structure to function?
2. can we label proteins to be used in supervised learning -> domains denoting functionality ?
3. Does similarity in structure imply similarity in function? i.e. If p1 and p2 are similar in structure, will they have same function?
