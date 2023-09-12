clear all
clc
%% main.m file to run the IEEEtripleCliped function
% Dicom image is given as input (As mostly medical images are in Dicom
% Format) with high dynamic range
image = dicomread('Subject_1.dcm');
%% Applying IEEE paper implementation code function
tcdhe_dwt_svd = IEEEtripleCliped(image);