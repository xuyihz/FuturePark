%% function
% MGT elevator
%
% Xu Yi, 24th April 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_elevator(fileID, iNO, iEL, CoC_elevator, Deg_elevator, facade_ele_R, levelZaxis_f, levelPstart1, elevatorColu_num, elevatorR, ~, ~, ROOF, Arc_itvl)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis_f = length(levelZaxis_f(:));

XYcoor_i = zeros(elevatorColu_num,2);   % 8个点(含电梯中间不落柱的1个点和圆心点)内筒XoY坐标第1(X)、2(Y)列。
elevatorColu_o_num = elevatorColu_num+4; % 内筒8个点,外筒12个点
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
    elevatorXY_o_Rtemp = sqrt(facade_ele_R(j)^2/2); % 2/5/8/11 正方形四角
    elevatorXY_o = [elevatorXY_o_temp,elevatorXYtemp; elevatorXY_o_Rtemp,elevatorXY_o_Rtemp; elevatorXYtemp,elevatorXY_o_temp;...
        -elevatorXYtemp,elevatorXY_o_temp; -elevatorXY_o_Rtemp,elevatorXY_o_Rtemp; -elevatorXY_o_temp,elevatorXYtemp;...
        -elevatorXY_o_temp,-elevatorXYtemp; -elevatorXY_o_Rtemp,-elevatorXY_o_Rtemp; -elevatorXYtemp,-elevatorXY_o_temp;...
        elevatorXYtemp,-elevatorXY_o_temp; elevatorXY_o_Rtemp,-elevatorXY_o_Rtemp; elevatorXY_o_temp,-elevatorXYtemp];
    
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
    for j = 1:elevatorColu_o_num % 外部12个柱子
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis_f(i));
    end
end
% 以直代曲部分
iNO_main_end = iNO; % 主节点终点备份，即以直代曲节点起点备份。
XY_Deg_num = zeros(lengthlevelZaxis_f,lengthXYcoor_f); % 以直代曲的各曲线的分隔节点数
P_start = zeros(1,2); P_end = zeros(1,2);
for i = levelPstart1:lengthlevelZaxis_f
    for j = 1:lengthXYcoor_f % 幕墙
        P_start(:) = XYcoor_o3(i,j,:);
        if j == lengthXYcoor_f
            P_end(:) = XYcoor_o3(i,1,:);
        else
            P_end(:) = XYcoor_o3(i,j+1,:);
        end
        [iNO, Deg_num] = MGT_arc_FE(fileID, iNO, levelZaxis_f(i), CoC_elevator, P_start, P_end, Arc_itvl);
        XY_Deg_num(i,j) = Deg_num;
    end
end
% 以直代曲部分
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
    for j = 1:lengthXYcoor_i/2 % 内筒4个点
        for k = 1:3 % 内筒每个柱悬挑三个梁
            iEL = iEL+1;
            iN1 = iNO + j + lengthXYcoor_all*(i-1); % 内筒4个角点
            iN2 = iNO + k + 3*(j-1) + lengthXYcoor_i + lengthXYcoor_all*(i-1); % 外筒8个点
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
% 幕墙外环梁
iNO_arc = iNO_main_end; % 初始化
for i = levelPstart1:lengthlevelZaxis_f	%
    for j = 1:lengthXYcoor_f
        iN1_bkp = iNO+j+lengthXYcoor_i+lengthXYcoor_all*(i-1); %
        if j ~= lengthXYcoor_f
            iN2_bkp = iN1_bkp+1;
        else % j = lengthXYcoor_f 时， 连接的是本环的第一个点，而不是上层内环的第一个点。
            iN2_bkp = iN1_bkp+1-lengthXYcoor_f;
        end
        
        if XY_Deg_num(i,j) == 1 % 即此处未进行以直代曲分隔
            iEL = iEL+1;
            iN1 = iN1_bkp;
            iN2 = iN2_bkp;
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 梁单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        else
            for k = 1:XY_Deg_num(i,j)
                iEL = iEL+1;
                if k == 1
                    iNO_arc = iNO_arc+1;
                    iN1 = iN1_bkp;
                    iN2 = iNO_arc;
                elseif k == XY_Deg_num(i,j)
                    iN1 = iNO_arc;
                    iN2 = iN2_bkp;
                else
                    iN1 = iNO_arc;
                    iNO_arc = iNO_arc+1;
                    iN2 = iNO_arc;
                end
                fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                    iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                    iN1, iN2,...    % 梁单元的两个节点号
                    ELE_ANGLE, ELE_iSUB);
            end
        end
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