%% Generate MGT file
% main M file
% 
% Xu Yi, 19th March 2018
% Xu Yi, 21st March 2018, revised

%%
close all; clear; clc;

%% 
fileID = fopen('FuturePark.mgt','w');   % Open or create new file for writing. Discard existing contents, if any.
addpath(genpath('coor_fun'))    % 搜索路径中加入coor_fun文件夹及其下所有文件夹
addpath(genpath('nurbs_fun'))

%% append initial conditions
MGT_init(fileID);

%% append model file
MGT_model(fileID);

%%
fprintf(fileID,'*ENDDATA');

%%
fclose('all');
