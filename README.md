# KPPM_KneePathologyPainMapping
Matlab code for the evaluation of Voxel-based Knee Pathology-Pain Mapping (KPPM) as described in the publication 'Arthofer et al. (2020) An anatomical atlas of the knee and voxel-based knee pathology-pain mapping'. Please cite if you use the code.

We implemented two methods: permutation-based thresholding with binary logistic regression and Mann-Whitney testing. Lesion and pain data can be simulated or provided as input.

Using the code:
1.  Specify:
- whether data is simulated or your own data provided
- the type of test
- an image to provide as reference space (for example this could be an MNI template used for spatial normalisation or your group template generated from your own data)
- an output folder
- for BLR also specify the number of permutations
    
    
