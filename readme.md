# Documentation for Xplore Autumn 2021 project


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