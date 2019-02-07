clc;
clear;

addpath(genpath('path-to-{CollaborativeBrainDecomposition}'));


candidateLstFile = 'path-to-output/init_cand_lst.txt';
outDir = 'path-to-output/robustInit';
K = 17;

initV = selRobustInit(candidateLstFile,K,outDir);
