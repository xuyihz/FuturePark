%% function
% MGT tower 3 facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_facade_S3(fileID, iNO, iEL, column_num, CoC_tower, Deg_tower, towerS_column_coor, facade_tower_R, levelZaxis, levelPstart, ~, ~)    %注意这里的levelPstart是1x2数组
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

levelPstart1 = levelPstart(1);
lengthlevelZaxis = length(levelZaxis(:));
xtemp = towerS_column_coor(1); ytemp = towerS_column_coor(2);
XYcoor_o3 = zeros(lengthlevelZaxis,column_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。注意这里是三维数组。(Z方向幕墙有变化)

Rec_R = sqrt( xtemp^2 + ytemp^2 ); % 4个柱子形成的矩形对角线的一半
for j = levelPstart1:lengthlevelZaxis
    Rtemp = facade_tower_R(j);
    XY_o_2x = xtemp / Rec_R * Rtemp;	% 幕墙上2点 x。
    XY_o_2y = ytemp / Rec_R * Rtemp;  % 幕墙上2点 y。
    
    XYcoor_o3(j,:,:) = [Rtemp, 0; XY_o_2x, XY_o_2y; 0, Rtemp; -XY_o_2x, XY_o_2y;...
        -Rtemp, 0; -XY_o_2x, -XY_o_2y; 0, -Rtemp; XY_o_2x, -XY_o_2y;];
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