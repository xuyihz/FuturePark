%% Generate MGTmisc file
% main M file
% 
% Xu Yi, 2018

%%
close all; clear; clc;

%% 
fileID = fopen('FuturePark_Misc.mgt','w');   % Open or create new file for writing. Discard existing contents, if any.
addpath(genpath('./coor_fun'))    % ����·���м�����һ��Ŀ¼��coor_fun�ļ��м����������ļ���

%% append model file
MGTmisc_model(fileID);

%%
fprintf(fileID,'*ENDDATA');

%%
fclose('all');
