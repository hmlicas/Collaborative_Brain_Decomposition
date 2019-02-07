clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

sbjListFile = 'path-to-data/fs_sbjLst.txt';
surfL = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/surf/lh.pial';
surfR = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/surf/rh.pial';
surfML = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/label/lh.Medial_wall.label';
surfMR = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/label/rh.Medial_wall.label';
prepDataFile = 'path-to-output-fs/prepData.mat';
outDir = 'path-to-output-fs/init';
spaR = 1;
vxI = 0;
ard = 0;
iterNum = 1000;
K = 17;
tNum = 120;
alpha = 2;
beta = 10;
resId = 'fmri_fs';

deployFuncInit_surf_fs(sbjListFile,surfL,surfR,surfML,surfMR,prepDataFile,outDir,spaR,vxI,ard,iterNum,K,tNum,alpha,beta,resId);
