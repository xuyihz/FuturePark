%% function
% MGT stairs
%
% Xu Yi, 25th March 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_stair(fileID, iNO, iEL, CoC_stair, Deg_stair, levelZaxis, levelPstart1, stairColu_num, stairL, stairW, stairB, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;

XYcoor_i = zeros(stairColu_num,2);   % 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o = zeros(stairColu_num*2,2); % 外表皮XoY坐标第1(X)、2(Y)列。

stairXY = [stairL/2, stairW/2; -stairL/2, stairW/2; -stairL/2, -stairW/2; stairL/2, -stairW/2]; % 原始坐标，未旋转，未转到整体坐标系
for i = 1:stairColu_num   % 尝试向量化
    [XYcoor_i(i,:)] = coorTrans(stairXY(i,:), Deg_stair); % 内筒
end
% 外筒待定stairXY2 暂定从内筒外伸stairB.
stairXY2 = [stairL/2+stairB, stairW/2; stairL/2, stairW/2+stairB; -stairL/2, stairW/2+stairB; -stairL/2-stairB, stairW/2;...
    -stairL/2-stairB, -stairW/2; -stairL/2, -stairW/2-stairB; stairL/2, -stairW/2-stairB; stairL/2+stairB, -stairW/2];
for i = 1:stairColu_num*2   % 尝试向量化
    [XYcoor_o(i,:)] = coorTrans(stairXY2(i,:), Deg_stair); % 外筒
end
% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_stair(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_stair(2);
XYcoor_o(:,1) = XYcoor_o(:,1) + CoC_stair(1);
XYcoor_o(:,2) = XYcoor_o(:,2) + CoC_stair(2);
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for k = 1:2 % 内筒，外筒
        if k == 1 % 节点编号规则：从0度角开始逆时针；先每层内筒，再每层外筒；从下到上。
            for j = 1:stairColu_num % 内部4个柱子
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
            end
        else
            for j = 1:stairColu_num*2 % 外部8个柱子
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_o(j,1),XYcoor_o(j,2),levelZaxis(i));
            end
        end
    end
end
lengthXYcoor2 = stairColu_num+stairColu_num*2;  % 每层的节点数，其中内部4个点，外部8个点。
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
    for j = 1:stairColu_num	% 每层内筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor2*(i-1);
        iN2 = iN1+lengthXYcoor2;
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
for i = levelPstart1:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数 % levelPstart 第几层开始停车，即下几层开敞
    for j = 1:stairColu_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+(stairColu_num+j)+lengthXYcoor2*(i-1); % 此行与内筒不同，多了 +stairN_num/2
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
fprintf(fileID,'; 楼梯长向主梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis-1)	% 由于有斜段，故这里要-1
    if rem(i,2) ~= 0    % 奇数层 % 控制点，即两个斜段起点
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor2+1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % 偶数层
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor2-1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+stairColu_num*2-1;
    end
    if i < levelPstart1 % 考虑底层无幕墙
        k_end = 2;
    else
        k_end = 4;
    end
    for k = 1:k_end % 梁有两根各两段 % 斜段+平段
        if k == 1
            iN1 = iNcon1; iN2 = iNcon3;
        elseif k == 2
            iN1 = iNcon2; iN2 = iNcon4;
        elseif k == 3
            iN1 = iNcon3; iN2 = iNcon5;
        elseif k == 4
            iN1 = iNcon4; iN2 = iNcon6;
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
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
    else % 偶数层
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+1;
    end
    if i < levelPstart1 % 考虑底层无幕墙
        k_end = 1;
    else
        k_end = 3;
    end
    for k = 1:k_end % 梁有三段
        if k == 1
            iN1 = iNcon1;
            iN2 = iNcon2;
        elseif k == 2
            iN1 = iNcon1;
            iN2 = iN1 + 5;
        elseif k == 3
            iN1 = iNcon2;
            iN2 = iN1 + 7;
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
for i = levelPstart1:(lengthlevelZaxis-1)	% 由于有斜段，故这里要-1
    if rem(i,2) ~= 0    % 奇数层 % 控制点，即两个斜段起点
        iNcon1 = iNO+6+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+5;
        iNcon3 = iNcon1+lengthXYcoor2+1;
        iNcon4 = iNcon2+lengthXYcoor2-1;
        iNcon5 = iNcon3+1;
        iNcon6 = iNcon4-1;
    else % 偶数层
        iNcon1 = iNO+7+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor2-1;
        iNcon4 = iNcon2+lengthXYcoor2+1;
        iNcon5 = iNcon3-1;
        iNcon6 = iNcon4+1;
    end
    for k = 1:5 % 梁有长向两根各两段，宽向一根一段 % 斜段+平段
        if k == 1
            iN1 = iNcon1; iN2 = iNcon3;
        elseif k == 2
            iN1 = iNcon2; iN2 = iNcon4;
        elseif k == 3
            iN1 = iNcon3; iN2 = iNcon5;
        elseif k == 4
            iN1 = iNcon4; iN2 = iNcon6;
        elseif k == 5
            iN1 = iNcon5; iN2 = iNcon6;
        end
        iEL = iEL+1;
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
for i = 1:(lengthlevelZaxis-1) % 由于有斜段，故这里要-1
    if rem(i,2) ~= 0    % 奇数层 % 控制点，即两个斜段起点
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor2+1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % 偶数层
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor2-1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+stairColu_num*2-1;
    end
    if i < levelPstart1 % 考虑底层无幕墙
        k_end = 1;
    else
        k_end = 2;
    end
    for k = 1:k_end % 两块板 % 斜段+平段
        if k == 1
            iN1 = iNcon1; iN2 = iNcon3; iN3 = iNcon4; iN4 = iNcon2;
        elseif k == 2
            iN1 = iNcon3; iN2 = iNcon5; iN3 = iNcon6; iN4 = iNcon4;
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2, iN3, iN4,...    % 板单元的四个节点号
            ELE_iSUB, ELE_iWID);
    end
end
fprintf(fileID,'\n');

%% FLOORLOAD
fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');

LTNAME = ROOF; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0; % 楼梯荷载 目前采用了屋面荷载4/3.5，板厚0.因每2.1一块板，故可能和4.2一层7/3.5差不多。(待复核)
DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';

iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis-1) % 由于有斜段，故这里要-1
    if rem(i,2) ~= 0    % 奇数层 % 控制点，即两个斜段起点
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor2+1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % 偶数层
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor2-1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+stairColu_num*2-1;
    end
    if i < levelPstart1 % 考虑底层无幕墙
        k_end = 1;
    else
        k_end = 2;
    end
    for k = 1:k_end % 两块板 % 斜段+平段
        if k == 1
            iN1 = iNcon1; iN2 = iNcon3; iN3 = iNcon4; iN4 = iNcon2;
        elseif k == 2
            iN1 = iNcon3; iN2 = iNcon5; iN3 = iNcon6; iN4 = iNcon4;
        end
        fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
    end
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % 初始化iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+stairColu_num);
CONSTRAINT = 111111; % 6个自由度全约束
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end