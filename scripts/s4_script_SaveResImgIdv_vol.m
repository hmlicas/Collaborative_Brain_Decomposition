clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

resFileName = 'path-to-output-vol/res/final_UV.mat';
maskName = 'path-to-data/mask.nii.gz';
outDir = 'path-to-output-vol/fig_idv';
saveFig = 1;
refNiiName = 'path-to-data/mask.nii.gz';  % reference image
sbjId = 1;

func_saveVolRes2Nii_idv(resFileName,maskName,sbjId,outDir,saveFig,refNiiName);
