clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));

% for volumetric data
maskFile = 'path-to-mask-file/mask.nii.gz';
maskNii = load_untouch_nii(maskFile);

gNb = createPrepData('volumetric', maskNii.img, 1);

% for surface data
surfL = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/surf/lh.pial';
surfR = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/surf/rh.pial';
surfML = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/label/lh.Medial_wall.label';
surfMR = 'path-to-{CollaborativeBrainDecomposition}/lib/freesurfer/subjects/fsaverage5/label/rh.Medial_wall.label';

[surfStru, surfMask] = getFsSurf(surfL, surfR, surfML, surfMR);

gNb = createPrepData('surface', surfStru, 1, surfMask);

% save gNb into file for later use
prepDataName = 'path-to-output/test_CreatePrepData.mat';
save(prepDataName, 'gNb');
