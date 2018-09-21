%% Generate MGTmisc file
% main M file
% 
% Xu Yi, 2018

%%
close all; clear; clc;

%% 
fileID = fopen('FuturePark_Misc.mgt','w');   % Open or create new file for writing. Discard existing contents, if any.
addpath(genpath('./coor_fun'))    % 搜索路径中加入上一级目录中coor_fun文件夹及其下所有文件夹

%% append model file
MGTmisc_model(fileID);

%%
fprintf(fileID,'*ENDDATA');

%%
fclose('all');
