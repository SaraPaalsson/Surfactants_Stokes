% close all; 
clear all
clc

fileName = '../output/data/swissroll3F_outdata_t0p1.mat';
load(fileName)


z = simdata.drops.z;
idx = simdata.drops.idx;

findClosestDist(z,idx)