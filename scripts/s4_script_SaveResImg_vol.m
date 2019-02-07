clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

resFileName = 'path-to-output-vol/res/res_cen.mat';
maskName = 'path-to-data/mask.nii.gz';
outDir = 'path-to-output-vol/fig';
saveFig = 1;
refNiiName = 'path-to-data/mask.nii.gz';  % reference image

func_saveVolRes2Nii(resFileName,maskName,outDir,saveFig,refNiiName);
