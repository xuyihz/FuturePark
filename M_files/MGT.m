%% Generate MGT file
% main M file
% 
% Xu Yi, 19th March 2018
% Xu Yi, 21st March 2018, revised

%%
close all; clear; clc;

%% 
fileID = fopen('FuturePark.mgt','w');   % Open or create new file for writing. Discard existing contents, if any.

%% append initial conditions
MGT_init(fileID);

%% append material conditions
MGT_mat(fileID);	% 目前只定义了Q345, C30

%% append model file
MGT_model(fileID);

%%
fprintf(fileID,'*ENDDATA');

%%
fclose('all');
