%% function
% MGT stairs
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_stair(fileID, iNO, iEL, CoC_stair, Deg_stair, facade_stair_R, levelZaxis_f, levelPstart1, stairColu_num, stairL, stairW, ~, ~, ~, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

lengthlevelZaxis_f = length(levelZaxis_f(:));

XYcoor_i = zeros(stairColu_num,2);   % 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o3 = zeros(lengthlevelZaxis_f,stairColu_num*3,2);	% 这里*3，外筒XoY坐标第1(X)、2(Y)列。注意这里是三维数组。(Z方向幕墙有变化)

% 内筒点
stairXYtemp1 = [stairL/2, stairW/2; -stairL/2, stairW/2]; % 原始坐标，未旋转，未转到整体坐标系
stairXYtemp1Deg = coorDeg([0,0], [1,0], stairXYtemp1(1,:)); % 内筒第一点与X轴的角度
stairXYtemp2 = zeros(length(stairXYtemp1),2);
for i = 1:length(stairXYtemp1)
    stairXYtemp2(i,:) = stairXYtemp1((length(stairXYtemp1)-i+1),:); % 逆时针编号
end
for i = 1:length(stairXYtemp1)
    stairXYtemp2(i,2) = -stairXYtemp1((length(stairXYtemp1)-i+1),2); % 沿X轴对称
end
stairXY_i = [stairXYtemp1; stairXYtemp2];

for i = 1:stairColu_num   % 尝试向量化
    [XYcoor_i(i,:)] = coorTrans(stairXY_i(i,:), Deg_stair); % 内筒
end
% 幕墙外筒
for j = levelPstart1:lengthlevelZaxis_f
    XYcoor_o_x = sqrt(facade_stair_R(j)^2 - (stairW/2)^2);  % 幕墙外圈点当y=±stairW/2时，x的绝对值
    XYcoor_o_y = sqrt(facade_stair_R(j)^2 - (stairL/2)^2);  % 幕墙外圈点当x=±stairL/2时，y的绝对值
    
    stairXYtemp = coorTrans([facade_stair_R(j), 0], -stairXYtemp1Deg); % 斜梁端点坐标。
    
    stairXY_o_temp1 = [XYcoor_o_x, stairW/2; stairXYtemp; stairL/2, XYcoor_o_y];  % 3个点的模块
    stairXY_o_temp2 = zeros(length(stairXY_o_temp1),2);
    for i = 1:length(stairXY_o_temp1)
        stairXY_o_temp2(i,:) = stairXY_o_temp1((length(stairXY_o_temp1)-i+1),:); % 逆时针编号
    end
    for i = 1:length(stairXY_o_temp1)
        stairXY_o_temp2(i,1) = -stairXY_o_temp1((length(stairXY_o_temp1)-i+1),1); % 沿Y轴对称
    end
    stairXY_o_temp12 = [stairXY_o_temp1; stairXY_o_temp2];
    stairXY_o_temp4 = zeros(length(stairXY_o_temp12),2);
    for i = 1:length(stairXY_o_temp12)
        stairXY_o_temp4(i,:) = stairXY_o_temp12((length(stairXY_o_temp12)-i+1),:); % 逆时针编号
    end
    for i = 1:length(stairXY_o_temp12)
        stairXY_o_temp4(i,2) = -stairXY_o_temp12((length(stairXY_o_temp12)-i+1),2); % 沿X轴对称
    end
    stairXY_o = [stairXY_o_temp1; stairXY_o_temp2; stairXY_o_temp4]; % 幕墙外筒各点坐标
    
    for i = 1:length(stairXY_o)   % 尝试向量化
        [XYcoor_o3(j,i,:)] = coorTrans(stairXY_o(i,:), Deg_stair); % 幕墙外筒点坐标 % 已旋转整体角度
    end
end

% 局部坐标系 转换至 整体坐标系
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_stair(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_stair(2);

for i = 1:lengthlevelZaxis_f  % length(A(:)) A向量元素个数
    % 内筒，外筒 % 节点编号规则：从0度角开始逆时针；先每层内筒，再每层外筒；从下到上。
    for j = 1:stairColu_num*3 % 外部12个柱子
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis_f(i));
    end
end
iNO_end = iNO;
iEL_end = iEL;
fprintf(fileID,'\n');
end