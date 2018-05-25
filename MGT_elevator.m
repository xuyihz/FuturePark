%% function
% MGT elevator
%
% Xu Yi, 24th April 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_elevator(fileID, iNO, iEL, CoC_elevator, Deg_elevator, facade_ele_R, levelZaxis_f, levelPstart1, elevatorColu_num, elevatorR, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis_f = length(levelZaxis_f(:));

XYcoor_i = zeros(elevatorColu_num,2);   % 8个点(含电梯中间不落柱的1个点和圆心点)内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o3 = zeros(lengthlevelZaxis_f,elevatorColu_num,2);	% 8个点，外筒XoY坐标第1(X)、2(Y)列。注意这里是三维数组。(Z方向幕墙有变化)
elevatorColu_o_num = elevatorColu_num; % 内外筒都8个点

% 内筒点
elevatorXYtemp = elevatorR/sqrt(2);
elevatorXY = [elevatorXYtemp,elevatorXYtemp; -elevatorXYtemp,elevatorXYtemp;...
    -elevatorXYtemp,-elevatorXYtemp; elevatorXYtemp,-elevatorXYtemp;...
    elevatorXYtemp,0; -elevatorXYtemp,0; 0,elevatorXYtemp; 0,0]; % 原始坐标，未旋转，未转到整体坐标系
for i = 1:elevatorColu_num   % 尝试向量化
    XYcoor_i(i,:) = coorTrans(elevatorXY(i,:), Deg_elevator); % 内筒
end
% 外筒elevatorXY2 % 外筒需要根据幕墙外表皮曲线定位 % 需有沿高度的循环(或向量化)
for j = levelPstart1:lengthlevelZaxis_f
    elevatorXY_o_temp = sqrt(facade_ele_R(j)^2 - elevatorXYtemp^2);
    elevatorXY_o = [elevatorXY_o_temp,elevatorXYtemp; elevatorXYtemp,elevatorXY_o_temp;...
        -elevatorXYtemp,elevatorXY_o_temp; -elevatorXY_o_temp,elevatorXYtemp;...
        -elevatorXY_o_temp,-elevatorXYtemp; -elevatorXYtemp,-elevatorXY_o_temp;...
        elevatorXYtemp,-elevatorXY_o_temp; elevatorXY_o_temp,-elevatorXYtemp];
    
    for i = 1:elevatorColu_o_num   % 尝试向量化
        [XYcoor_o3(j,i,:)] = coorTrans(elevatorXY_o(i,:), Deg_elevator); % 内筒
    end
end
% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_elevator(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_elevator(2);
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_elevator(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_elevator(2);

lengthXYcoor_i = length(XYcoor_i); % 内筒每层节点数
lengthXYcoor_f = length(elevatorXY_o); % 幕墙每层节点数
lengthXYcoor_all = lengthXYcoor_i + lengthXYcoor_f;  % 每层节点数备份

for i = 1:lengthlevelZaxis_f  % length(A(:)) A向量元素个数
    % 内筒，外筒 % 节点编号规则：先每层内筒，再每层外筒；从下到上。
    for j = 1:elevatorColu_num % 内部7个柱子(有一个不落柱)
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis_f(i));
    end
    for j = 1:elevatorColu_o_num % 外部10个柱子
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis_f(i));
    end
end
iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_colu = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 内筒柱；iPRO = 2 截面编号2。
fprintf(fileID,'; 内筒柱\n');
ELE_iPRO = 2;
iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis_f-1)	% length(A(:)) A向量元素个数
    for j = 1:lengthXYcoor_i	% 每层内筒的节点数
        if j==7 || j==8 % 电梯中间节点不落柱
        else
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcoor_all*(i-1);
            iN2 = iN1+lengthXYcoor_all;
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 柱单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        end
    end
end

% 外筒柱；iPRO = 1 截面编号1。
fprintf(fileID,'; 外筒柱\n');
ELE_iPRO = 1;
iNO = iNO_init; % 初始化iNO
for i = levelPstart1:(lengthlevelZaxis_f-1)	% length(A(:)) A向量元素个数 % levelPstart 第几层开始停车，即下几层开敞
    for j = 1:lengthXYcoor_f	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+(lengthXYcoor_i+j)+lengthXYcoor_all*(i-1); % 此行与内筒不同，多了 +elevatorColu_num
        iN2 = iN1+lengthXYcoor_all;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 柱单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) beams 布置与停车筒不同，且均为同一截面
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 主梁；iPRO = 3 截面编号3。
fprintf(fileID,'; 电梯筒内筒主梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
iN_table = [1,2; 2,3; 3,4; 4,1; 5,6; 7,8]; % 表驱动 % 每行为梁两端节点号
for i = 1:lengthlevelZaxis_f	% 此行i起始为2.即二层开始有。
    for k = 1:length(iN_table) % 梁有6根
        iEL = iEL+1;
        iN1 = iNO + iN_table(k,1) + lengthXYcoor_all*(i-1);
        iN2 = iNO + iN_table(k,2) + lengthXYcoor_all*(i-1);
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

% 悬臂梁；iPRO = 3 截面编号3。
fprintf(fileID,'; 悬臂梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = levelPstart1:lengthlevelZaxis_f	%
    for j = 1:lengthXYcoor_f/2 % 内筒4个点
        for k = 1:2 % 内筒每个柱悬挑两个梁
            iEL = iEL+1;
            iN1 = iNO + j + lengthXYcoor_all*(i-1); % 内筒4个角点
            iN2 = iNO + k + 2*(j-1) + lengthXYcoor_i + lengthXYcoor_all*(i-1); % 外筒8个点
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 梁单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        end
    end
end

% 环次梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 电梯筒外环次梁\n');
ELE_iPRO = 4;
iNO = iNO_init; % 初始化iNO
for i = levelPstart1:lengthlevelZaxis_f	%
    for j = 1:lengthXYcoor_f
        iEL = iEL+1;
        iN1 = iNO + j + lengthXYcoor_i + lengthXYcoor_all*(i-1);
        if j == lengthXYcoor_f
            iN2 = iN1 + 1 - lengthXYcoor_f;
        else
            iN2 = iN1 + 1;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % 初始化iNO
% 1~6 因第7个节点没有立柱子，故这里跳过7
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+lengthXYcoor_i-2);
CONSTRAINT = 111111; % 6个自由度全约束
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end