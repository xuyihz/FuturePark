%% function
% MGT elevator
%
% Xu Yi, 24th April 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_elevator(fileID, iNO, iEL, CoC_elevator, Deg_elevator, facade_ele_R, levelZaxis, levelPstart1, elevatorColu_num, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
elevatorColu_o_num = elevatorColu_num+2;
lengthlevelZaxis = length(levelZaxis(:));
XYcoor_i = zeros(elevatorColu_num,2);   % 7个 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o3 = zeros(lengthlevelZaxis,elevatorColu_o_num,2); % 10个 外表皮XoY坐标第1(X)、2(Y)列。注意这里是三维数组。(Z方向幕墙有变化)

ele_width = 2500; ele_shift = 400; ele_depth = ele_width + ele_shift; % 内筒柱定位，待确定
str_length = 2200; str_width = 2800;
elevatorXY = [-ele_width, ele_depth; ele_width, ele_depth; -ele_width, ele_shift; ele_width, ele_shift;...
                -str_length, -str_width; str_length, -str_width; 0, ele_depth; 0, ele_shift]; % 原始坐标，未旋转，未转到整体坐标系
for i = 1:elevatorColu_num   % 尝试向量化
    XYcoor_i(i,:) = coorTrans(elevatorXY(i,:), Deg_elevator); % 内筒
end
% 外筒elevatorXY2 % 外筒需要根据幕墙外表皮曲线定位 % 需有沿高度的循环(或向量化)
for j = levelPstart1:lengthlevelZaxis
    ele2_depth = sqrt(abs(facade_ele_R(j)^2 - ele_width^2));
    ele2_width1 = sqrt(abs(facade_ele_R(j)^2 - ele_depth^2));
    ele2_width2 = sqrt(abs(facade_ele_R(j)^2 - ele_shift^2));
    ele2_width3 = sqrt(abs(facade_ele_R(j)^2 - str_width^2));
    [ele2_strX, ele2_strY] = coorLxCp([0,0], elevatorXY(3,:), elevatorXY(5,:), facade_ele_R(j));
    
    elevatorXY2 = [-ele_width, ele2_depth; ele_width, ele2_depth; -ele2_width1, ele_depth; ele2_width1, ele_depth;...
        -ele2_width2, ele_shift; ele2_width2, ele_shift; -ele2_width3, -str_width; ele2_width3, -str_width;...
        ele2_strX, ele2_strY; -ele2_strX, ele2_strY];
    for i = 1:elevatorColu_o_num   % 尝试向量化
        [XYcoor_o3(j,i,:)] = coorTrans(elevatorXY2(i,:), Deg_elevator); % 内筒
    end
end
% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_elevator(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_elevator(2);
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_elevator(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_elevator(2);
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    % 内筒，外筒 % 节点编号规则：先每层内筒，再每层外筒；从下到上。
    for j = 1:elevatorColu_num % 内部7个柱子
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
    end
    for j = 1:elevatorColu_o_num % 外部10个柱子
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));
    end
end
lengthXYcoor2 = elevatorColu_num+elevatorColu_o_num;  % 每层的节点数，其中内部47个点，外部10个点。
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
for i = 1:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数
    for j = 1:elevatorColu_num	% 每层内筒的节点数
        if j == 7   % 电梯中间那个节点不落柱
        else
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcoor2*(i-1);
            iN2 = iN1+lengthXYcoor2;
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
for i = levelPstart1:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数 % levelPstart 第几层开始停车，即下几层开敞
    for j = 1:elevatorColu_o_num	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+(elevatorColu_num+j)+lengthXYcoor2*(i-1); % 此行与内筒不同，多了 +elevatorColu_num
        iN2 = iN1+lengthXYcoor2;
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
fprintf(fileID,'; 电梯主梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = 2:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1; 此行i起始为2.即二层开始有。
    iNcon = iNO+lengthXYcoor2*(i-1);
    for k = 1:7 % 梁有7根
        switch k
            case 1
                iN1 = iNcon+1; iN2 = iNcon+3;
            case 2
                iN1 = iNcon+7; iN2 = iNcon+8;
            case 3
                iN1 = iNcon+2; iN2 = iNcon+4;
            case 4
                iN1 = iNcon+1; iN2 = iNcon+7;
            case 5
                iN1 = iNcon+7; iN2 = iNcon+2;
            case 6
                iN1 = iNcon+3; iN2 = iNcon+8;
            case 7
                iN1 = iNcon+8; iN2 = iNcon+4;
        end        
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

% 主梁；iPRO = 3 截面编号3。
fprintf(fileID,'; 楼梯长向主梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis-1)	% 由于有斜段，故这里要-1
    iNcon3 = iNO+3+lengthXYcoor2*(i-1);	% 节点3
    iNcon4 = iNcon3+1;                  % 节点4
    iNcon5 = iNcon3+2;                  % 节点5
    iNcon6 = iNcon3+3;                  % 节点6
    if rem(i,2) ~= 0    % 奇数层 % 控制点，即两个斜段起点
        iNcon4 = iNcon4+lengthXYcoor2; % 暂定节点3\5起点
        iNcon6 = iNcon6+lengthXYcoor2;
    else % 偶数层
        iNcon3 = iNcon3+lengthXYcoor2;
        iNcon5 = iNcon5+lengthXYcoor2;
    end
    for k = 1:2
        switch k
            case 1
                iN1 = iNcon3;
                iN2 = iNcon4;
            case 2
                iN1 = iNcon5;
                iN2 = iNcon6;
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'; 楼梯宽向主梁\n');
iNO = iNO_init; % 初始化iNO
for i = 1:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1 % 每层一根贯穿宽向主梁
    if rem(i,2) ~= 0    % 奇数层 % 控制点，即两个内筒悬挑起点
        iN1 = iNO+3+lengthXYcoor2*(i-1);
        iN2 = iN1+2;
    else % 偶数层
        iN1 = iNO+4+lengthXYcoor2*(i-1);
        iN2 = iN1+2;
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % 梁单元的两个节点号
        ELE_ANGLE, ELE_iSUB);
end

% 悬臂梁；iPRO = 3 截面编号3。
fprintf(fileID,'; 悬臂梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = levelPstart1:lengthlevelZaxis	%
    iNcon = iNO+lengthXYcoor2*(i-1);
    iNcon_o = iNcon + elevatorColu_num; % 外筒节点起点
    for k = 1:10 % 悬臂梁有10根
        switch k
            case 1
                iN1 = iNcon+1; iN2 = iNcon_o+1;
            case 2
                iN1 = iNcon+1; iN2 = iNcon_o+3;
            case 3
                iN1 = iNcon+2; iN2 = iNcon_o+2;
            case 4
                iN1 = iNcon+2; iN2 = iNcon_o+4;
            case 5
                iN1 = iNcon+3; iN2 = iNcon_o+5;
            case 6
                iN1 = iNcon+4; iN2 = iNcon_o+6;
            case 7
                iN1 = iNcon+5; iN2 = iNcon_o+7;
            case 8
                iN1 = iNcon+5; iN2 = iNcon_o+9;
            case 9
                iN1 = iNcon+6; iN2 = iNcon_o+8;
            case 10
                iN1 = iNcon+6; iN2 = iNcon_o+10;
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

% 环次梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 环形次梁\n');
ELE_iPRO = 4;
iNO = iNO_init; % 初始化iNO
% 外环梁 % 参考楼梯长向主梁
fprintf(fileID,';   外环梁\n');
for i = levelPstart1:lengthlevelZaxis	%
    iNcon_o = iNO+lengthXYcoor2*(i-1) + elevatorColu_num; % 外筒节点起点
    for k = 1:10 % 梁有10段
        switch k
            case 1
                iN1 = iNcon_o+1; iN2 = iNcon_o+2;
            case 2
                iN1 = iNcon_o+1; iN2 = iNcon_o+3;
            case 3
                iN1 = iNcon_o+2; iN2 = iNcon_o+4;
            case 4
                iN1 = iNcon_o+3; iN2 = iNcon_o+5;
            case 5
                iN1 = iNcon_o+4; iN2 = iNcon_o+6;
            case 6
                iN1 = iNcon_o+5; iN2 = iNcon_o+7;
            case 7
                iN1 = iNcon_o+6; iN2 = iNcon_o+8;
            case 8
                iN1 = iNcon_o+7; iN2 = iNcon_o+9;
            case 9
                iN1 = iNcon_o+8; iN2 = iNcon_o+10;
            case 10
                iN1 = iNcon_o+9; iN2 = iNcon_o+10;
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(planner) floor 由于该4个点不在一个平面，故无法施加楼面荷载
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_floor = iEL;
ELE_TYPE = 'PLATE'; ELE_iMAT = 2; ELE_iSUB = 2; ELE_iWID = 0; % iMAT = 2材料混凝土C30 % iSUB = 2 薄板

% 板厚1；iPRO = 2 截面编号2。
fprintf(fileID,'; 1厚板楼梯板\n');
ELE_iPRO = 2;
iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis-1) % 由于有斜段，故这里要-1
    if rem(i,2) ~= 0    % 奇数层 % 控制点，即两个斜段起点
        iN1 = iNO+3+lengthXYcoor2*(i-1);
        iN2 = iN1+2;
        iN3 = iN1+3+lengthXYcoor2;
        iN4 = iN1+1+lengthXYcoor2;
    else % 偶数层
        iN1 = iNO+6+lengthXYcoor2*(i-1);
        iN2 = iN1-2;
        iN3 = iN1-3+lengthXYcoor2;
        iN4 = iN1-1+lengthXYcoor2;
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2, iN3, iN4,...    % 板单元的四个节点号
        ELE_iSUB, ELE_iWID);
end
fprintf(fileID,'\n');

%% FLOORLOAD 由于该4个点不在一个平面，故无法施加楼面荷载
% fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
% fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');
% 
% LTNAME = ROOF; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0; % 楼梯荷载 目前采用了屋面荷载4/3.5，板厚0.因每2.1一块板，故可能和4.2一层7/3.5差不多。(待复核)
% DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';
% 
% iNO = iNO_init; % 初始化iNO
% for i = 1:(lengthlevelZaxis-1) % 由于有斜段，故这里要-1
%     if rem(i,2) ~= 0    % 奇数层 % 控制点，即两个斜段起点
%         iN1 = iNO+3+lengthXYcoor2*(i-1);
%         iN2 = iN1+2;
%         iN3 = iN1+3+lengthXYcoor2;
%         iN4 = iN1+1+lengthXYcoor2;
%     else % 偶数层
%         iN1 = iNO+6+lengthXYcoor2*(i-1);
%         iN2 = iN1-2;
%         iN3 = iN1-3+lengthXYcoor2;
%         iN4 = iN1-1+lengthXYcoor2;
%     end
%     fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
%         LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
%         iN1, iN2, iN3, iN4);
% end
% fprintf(fileID,'\n');

%%
iEL_end = iEL;

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % 初始化iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+elevatorColu_num);
CONSTRAINT = 111111; % 6个自由度全约束
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end