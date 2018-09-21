%% function
% MGT tower 2 facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_facade_S2(fileID, iNO, iEL, column_num, CoC_tower, Deg_tower, towerS_column_coor, facade_tower_R, levelZaxis, levelPstart, ~, ~)    %注意这里的levelPstart是1x2数组
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

levelPstart1 = levelPstart(1);
lengthlevelZaxis = length(levelZaxis(:));
xtemp = towerS_column_coor(1); ytemp = towerS_column_coor(2);
XYcoor_o3 = zeros(lengthlevelZaxis,column_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。注意这里是三维数组。(Z方向幕墙有变化)

for j = levelPstart1:lengthlevelZaxis
    XY_o_1x = sqrt(facade_tower_R(j)^2 - ytemp^2);	% 幕墙上1点 x。
    XY_o_2y = sqrt(facade_tower_R(j)^2 - xtemp^2);  % 幕墙上2点 y。
    
    XYcoor_o3(j,:,:) = [XY_o_1x, ytemp; xtemp, XY_o_2y; -xtemp, XY_o_2y; -XY_o_1x, ytemp;...
        -XY_o_1x, -ytemp; -xtemp, -XY_o_2y; xtemp, -XY_o_2y; XY_o_1x, -ytemp;];
    if Deg_tower ~= 0
        for i = 1:column_num*2
            XYcoor_o3(j,i,:) = coorTrans(XYcoor_o3(j,i,:), Deg_tower);
        end
    end
end

% 局部坐标系 转换至 整体坐标系
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_tower(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_tower(2);

lengthXYcoor_f = column_num*2;  % 幕墙每层节点数

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for j = 1:lengthXYcoor_f % 幕墙
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));   % 外筒 X & Y
    end
end
iNO_end = iNO;
iEL_end = iEL;
fprintf(fileID,'\n');
end