% RUN THIS FILE TO RUN THE COMPLETE PROCEDURE 
% Rest files are stored in folder

clear; clc; close all;
addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('   Stage 1: ocv_soc.m ');
ocv_soc;

fprintf('\n   Stage 2: parameters.m \n');
parameters;

fprintf('\n   Stage 3: test_udds.m \n');
test_udds;

fprintf('\n   Stage 4: test_US06.m \n');
test_US06;
