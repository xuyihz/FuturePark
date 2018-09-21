%% function
% MGT elevator
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_elevator(fileID, iNO, iEL, CoC_elevator, Deg_elevator, facade_ele_R, levelZaxis_f, levelPstart1, elevatorColu_num, elevatorR, ~, ~, ~, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

lengthlevelZaxis_f = length(levelZaxis_f(:));

XYcoor_i = zeros(elevatorColu_num,2);   % 8个点(含电梯中间不落柱的1个点和圆心点)内筒XoY坐标第1(X)、2(Y)列。
elevatorColu_o_num = elevatorColu_num; % 内筒8个点,外筒8个点
XYcoor_o3 = zeros(lengthlevelZaxis_f,elevatorColu_o_num,2);	% 12个点，外筒XoY坐标第1(X)、2(Y)列。注意这里是三维数组。(Z方向幕墙有变化)

% 内筒点
elevatorXYtemp = elevatorR/sqrt(2); % 正方形
elevatorXY = [elevatorXYtemp,elevatorXYtemp; -elevatorXYtemp,elevatorXYtemp;...
    -elevatorXYtemp,-elevatorXYtemp; elevatorXYtemp,-elevatorXYtemp;...
    elevatorXYtemp,350; -elevatorXYtemp,350; 0,elevatorXYtemp; 0,350]; % 原始坐标，未旋转，未转到整体坐标系
for i = 1:elevatorColu_num   % 尝试向量化
    XYcoor_i(i,:) = coorTrans(elevatorXY(i,:), Deg_elevator); % 内筒
end
% 外筒elevatorXY2 % 外筒需要根据幕墙外表皮曲线定位 % 需有沿高度的循环(或向量化)
for j = levelPstart1:lengthlevelZaxis_f
    elevatorXY_o_temp = sqrt(facade_ele_R(j)^2 - elevatorXYtemp^2); %
    elevatorXY_o = [elevatorXY_o_temp,elevatorXYtemp; elevatorXYtemp,elevatorXY_o_temp;...
        -elevatorXYtemp,elevatorXY_o_temp; -elevatorXY_o_temp,elevatorXYtemp;...
        -elevatorXY_o_temp,-elevatorXYtemp; -elevatorXYtemp,-elevatorXY_o_temp;...
        elevatorXYtemp,-elevatorXY_o_temp; elevatorXY_o_temp,-elevatorXYtemp];
    
    for i = 1:elevatorColu_o_num   % 尝试向量化
        [XYcoor_o3(j,i,:)] = coorTrans(elevatorXY_o(i,:), Deg_elevator); % 内筒
    end
end
% 局部坐标系 转换至 整体坐标系
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_elevator(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_elevator(2);

for i = 1:lengthlevelZaxis_f  % length(A(:)) A向量元素个数
    % 内筒，外筒 % 节点编号规则：先每层内筒，再每层外筒；从下到上。
    for j = 1:elevatorColu_o_num % 外部8个柱子
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis_f(i));
    end
end
iNO_end = iNO;
iEL_end = iEL;
fprintf(fileID,'\n');
end