clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

sbjListFile = 'path-to-data/fs_sbjLst.txt';
medialWallFileL = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/label/lh.Medial_wall.label';
medialWallFileR = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/label/rh.Medial_wall.label';
prepDataFile = 'path-to-output-fs/prepData.mat';
outDir = 'path-to-output-fs/res';

resId = 'fmri_fs';
initName = 'path-to-output/robustInit/init.mat';
K = 17;
alphaS21 = 2;
alphaL = 10;
vxI = 0;
spaR = 1;
ard = 1;
eta = 1;
iterNum = 30;
calcGrp = 1;
parforOn = 0;

deployFuncMvnmfL21p1_func_surf_fs(sbjListFile,medialWallFileL,medialWallFileR,prepDataFile,outDir,resId,initName,K,alphaS21,alphaL,vxI,spaR,ard,eta,iterNum,calcGrp,parforOn);
