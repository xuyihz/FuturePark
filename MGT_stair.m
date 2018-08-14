%% function
% MGT stairs
%
% Xu Yi, 25th March 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_stair(fileID, iNO, iEL, CoC_stair, Deg_stair, facade_stair_R, levelZaxis_f, levelPstart1, stairColu_num, stairL, stairW, ~, ~, ROOF, Arc_itvl)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
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
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_stair(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_stair(2);
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_stair(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_stair(2);

lengthXYcoor_i = length(XYcoor_i); % 内筒每层节点数
lengthXYcoor_f = length(stairXY_o); % 幕墙每层节点数
lengthXYcoor_all = lengthXYcoor_i + lengthXYcoor_f;  % 每层节点数备份

for i = 1:lengthlevelZaxis_f  % length(A(:)) A向量元素个数
    % 内筒，外筒 % 节点编号规则：从0度角开始逆时针；先每层内筒，再每层外筒；从下到上。
    for j = 1:lengthXYcoor_i % 内部4个柱子
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis_f(i));
    end
    for j = 1:lengthXYcoor_f % 外部12个柱子
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
        [iNO, Deg_num] = MGT_arc_FE(fileID, iNO, levelZaxis_f(i), CoC_stair, P_start, P_end, Arc_itvl);
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
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_all*(i-1);
        iN2 = iN1+lengthXYcoor_all;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 柱单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

% 外筒柱；iPRO = 1 截面编号1。
fprintf(fileID,'; 幕墙柱\n');
ELE_iPRO = 1;
iNO = iNO_init; % 初始化iNO
for i = levelPstart1:(lengthlevelZaxis_f-1)	% length(A(:)) A向量元素个数 % levelPstart 第几层开始停车，即下几层开敞
    for j = 1:lengthXYcoor_f	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+(lengthXYcoor_i+j)+lengthXYcoor_all*(i-1); % 此行与内筒不同，多了 +stairN_num/2
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
fprintf(fileID,'; 楼梯主梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = 1:lengthlevelZaxis_f	%
    for j = 1:lengthXYcoor_i	% 每层内筒的节点数 其中2~3(建休息平台时建立platform.m)、4~1节点不连(休息平台净高要求)
        if j == 2 || j == 4
            % 跳过 2~3、4~1节点
        else
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcoor_all*(i-1);
            iN2 = iN1+1;    %
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 梁单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        end
    end
end

fprintf(fileID,'; 楼梯悬挑梁\n');
iNO = iNO_init; % 初始化iNO
for i = levelPstart1:lengthlevelZaxis_f % 楼梯都是从外筒伸出
    for j = 1:lengthXYcoor_i % 每层内筒的节点数
        for k = 1:(lengthXYcoor_f/lengthXYcoor_i) % 每层内筒对应的外筒节点数
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcoor_all*(i-1); % 此行为定位梁在塔楼的节点(内筒)
            iN2 = iNO+lengthXYcoor_i+(lengthXYcoor_f/lengthXYcoor_i)*(j-1)+k+lengthXYcoor_all*(i-1);    % 归到幕墙外筒第0点后，再定位到具体点
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 梁单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        end
    end
end

% 环形次梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 环形次梁\n');
ELE_iPRO = 4;
iNO = iNO_init; % 初始化iNO
% 外环梁
fprintf(fileID,';   幕墙外环梁\n');
iNO_arc = iNO_main_end; % 初始化
for i = levelPstart1:lengthlevelZaxis_f	% 此行与柱单元不同，柱单元为i-1;
    for j = 1:lengthXYcoor_f	% 每层外筒的节点数
        iN1_bkp = iNO+lengthXYcoor_i+j+lengthXYcoor_all*(i-1); %
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
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+lengthXYcoor_i);
CONSTRAINT = 111111; % 6个自由度全约束
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end