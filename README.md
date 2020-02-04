# KPPM_KneePathologyPainMapping
Matlab script for the evaluation of Voxel-based Knee Pathology-Pain Mapping (KPPM) as described in 'Arthofer et al. (2020) An anatomical atlas of the knee and voxel-based knee pathology-pain mapping'. Please cite if you use the code.

We implemented two methods: permutation-based thresholding with binary logistic regression and Mann-Whitney testing (currently not used). Lesion and pain data can be simulated or provided as input with a design file. Tested with Matlab R2018b.

## Using the code

Specify:
- whether data is simulated or your own data provided
- the type of test
- an image to provide the reference space (for example this could be an MNI template that was used for spatial normalisation or the group template generated in our work)
- an output folder
- for BLR also specify the number of permutations
- the design file in xlsx format
    
