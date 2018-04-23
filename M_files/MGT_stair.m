%% function
% MGT stairs
%
% Xu Yi, 25th March 2018
% Xu Yi, 28th March 2018, revised

%%
function [iNO_end, iEL_end] = MGT_stair(fileID, iNO, iEL, CoC_tower, levelZaxis, levelPstart, stairN_num, stairL, stairW, CAR, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

stairN_num = 8;
iNO_init = iNO;
XYcor = zeros(stairN_num,4);   % 内筒XoY坐标第1(X)、2(Y)列(1~4行)(注：5~8行均为0，不使用)，外表皮XoY坐标第3(X)、4(Y)列(1~8行)。

stairL = 3150; % 楼梯长，即沿踏步前进方向长
stairW = 3250; % 楼梯宽
stairXY = [stairL/2, stairW/2; -stairL/2, stairW/2; -stairL/2, -stairW/2; stairL/2, -stairW/2]; % 原始坐标，未旋转，未转到整体坐标系
stairDeg = -pi/6; % 楼梯旋转角度，暂定 顺时针转30(pi/6)
% 坐标系转换 function coorTrans
% x' = xcos(a) + ysin(a);
% y' = ycos(a) - xsin(a);
for i = 1:stairN_num/2   % 尝试向量化
    [XYcor(i,1), XYcor(i,2)] = coorTrans(stairXY(i,1), stairXY(i,2), stairDeg); % 内筒
end
% 外筒待定stairXY2 暂定从内筒外伸2000.
m2 = 2000;
stairXY2 = [stairL/2+m2, stairW/2; stairL/2, stairW/2+m2; -stairL/2, stairW/2+m2; -stairL/2-m2, stairW/2;...
    -stairL/2-m2, -stairW/2; -stairL/2, -stairW/2-m2; stairL/2, -stairW/2-m2; stairL/2+m2, -stairW/2];
for i = 1:stairN_num   % 尝试向量化
    [XYcor(i,3), XYcor(i,4)] = coorTrans(stairXY2(i,1), stairXY2(i,2), stairDeg); % 内筒
end
% 局部坐标系 转换至 整体坐标系
XYcor(:,1) = XYcor(:,1) + CoC_tower(1);
XYcor(:,2) = XYcor(:,2) + CoC_tower(2);
XYcor(:,3) = XYcor(:,3) + CoC_tower(1);
XYcor(:,4) = XYcor(:,4) + CoC_tower(2);
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for k = 1:2 % 内筒，外筒
        if k == 1 % 节点编号规则：从0度角开始逆时针；先每层内筒，再每层外筒；从下到上。
            for j = 1:stairN_num/2 % 内部4个柱子
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcor(j,1),XYcor(j,2),levelZaxis(i));
            end
        else
            for j = 1:stairN_num % 外部8个柱子
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcor(j,3),XYcor(j,4),levelZaxis(i));
            end
        end
    end
end
lengthXYcor2 = stairN_num/2+stairN_num;  % 每层的节点数，其中内部4个点，外部8个点。
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
    for j = 1:stairN_num/2	% 每层内筒的节点数 stairN_num/2
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcor2*(i-1);
        iN2 = iN1+lengthXYcor2;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 柱单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

% 外筒柱；iPRO = 1 截面编号1。
fprintf(fileID,'; 外筒柱\n');
ELE_iPRO = 1;
iNO = iNO_init; % 初始化iNO
for i = levelPstart:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数 % levelPstart 第几层开始停车，即下几层开敞
    for j = 1:stairN_num	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+(stairN_num/2+j)+lengthXYcor2*(i-1); % 此行与内筒不同，多了 +stairN_num/2
        iN2 = iN1+lengthXYcor2;
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

% 横向主梁；iPRO = 3 截面编号3。
fprintf(fileID,'; 横向主梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = levelPstart:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1
    for j = 1:stairN_num/2	% 每层内筒的节点数
        iN1 = iNO+j+lengthXYcor2*(i-1);
        for k = 1:2
            iN2 = iN1+stairN_num/2-(j-1)+(j-1)*2+(k-1); % 与停车筒的区别在这一行
            iEL = iEL+1;
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
% 内环梁
fprintf(fileID,';   内环梁\n');
for i = 2:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1; 此行与主梁不同，i起始为2.即二层开始有。
    for j = 1:stairN_num/2	% 每层内筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcor2*(i-1);
        if j ~= stairN_num/2
            iN2 = iN1+1;
        else % j = car_num 时， 连接的是本环的第一个点，而不是外环的第一个点。
            iN2 = iN1+1-stairN_num/2;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');
% 外环梁
fprintf(fileID,';   外环梁\n');
for i = levelPstart:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1;
    for j = 1:stairN_num	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+stairN_num/2+j+lengthXYcor2*(i-1); % 此行与内环梁不同，多加了stairN_num/2
        if j ~= stairN_num
            iN2 = iN1+1;
        else % j = car_num 时， 连接的是本环的第一个点，而不是上层内环的第一个点。
            iN2 = iN1+1-stairN_num;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(planner) floor
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_floor = iEL;
ELE_TYPE = 'PLATE'; ELE_iMAT = 2; ELE_iSUB = 2; ELE_iWID = 0; % iMAT = 2材料混凝土C30 % iSUB = 2 薄板

% 板厚1；iPRO = 2 截面编号2。
fprintf(fileID,'; 1厚板楼梯板\n');
ELE_iPRO = 2;
iNO = iNO_init; % 初始化iNO
for i = 2:lengthlevelZaxis % 此行同外环梁 % 此行与主梁不同，i起始为2.即二层开始有 % 因每层只定义了一块板，故没有j的循环
    iEL = iEL+1;
    iN1 = iNO+1+lengthXYcor2*(i-1); % 逆时针板四周四个点
    iN2 = iN1+1;
    iN3 = iN2+1;
    iN4 = iN3+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2, iN3, iN4,...    % 板单元的四个节点号
        ELE_iSUB, ELE_iWID);
end
fprintf(fileID,'\n');

%% FLOORLOAD
fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');

LTNAME = ROOF; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0; % 楼梯荷载 目前采用了屋面荷载4/3.5，板厚0.因每2.1一块板，故可能和4.2一层7/3.5差不多。(待复核)
DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';

iNO = iNO_init; % 初始化iNO
for i = 2:lengthlevelZaxis % 此行同外环梁 % 此行与主梁不同，i起始为2.即二层开始有 % 因每层只定义了一块板，故没有j的循环
    iN1 = iNO+1+lengthXYcor2*(i-1); % 逆时针板四周四个点
    iN2 = iN1+1;
    iN3 = iN2+1;
    iN4 = iN3+1;
    fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
        LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
        iN1, iN2, iN3, iN4);
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % 初始化iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+stairN_num/2);
CONSTRAINT = 111111; % 6个自由度全约束
fprintf(fileID,'   %s, %d, \n',...
                NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end