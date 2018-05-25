%% function
% MGT stairs
%
% Xu Yi, 25th March 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_stair_platform(fileID, iNO, iEL, CoC_stair, Deg_stair, levelZaxis, stairColu_num, stairL, stairW, stairB, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis = length(levelZaxis(:));

XYcoor_i = zeros(stairColu_num,2);   % 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_pf = XYcoor_i;	% 楼层平台/休息平台 XoY坐标第1(X)、2(Y)列。

% 内筒点
stairXYtemp1 = [stairL/2, stairW/2; -stairL/2, stairW/2]; % 原始坐标，未旋转，未转到整体坐标系
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
% 楼层平台/休息平台
stairXYtemp1 = [stairL/2+stairB, stairW/2; -stairL/2-stairB, stairW/2]; % 原始坐标，未旋转，未转到整体坐标系
stairXYtemp2 = zeros(length(stairXYtemp1),2);
for i = 1:length(stairXYtemp1)
    stairXYtemp2(i,:) = stairXYtemp1((length(stairXYtemp1)-i+1),:); % 逆时针编号
end
for i = 1:length(stairXYtemp1)
    stairXYtemp2(i,2) = -stairXYtemp1((length(stairXYtemp1)-i+1),2); % 沿X轴对称
end
stairXY_pf = [stairXYtemp1; stairXYtemp2];

for i = 1:stairColu_num   % 尝试向量化
    [XYcoor_pf(i,:)] = coorTrans(stairXY_pf(i,:), Deg_stair); % 内筒
end

% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_stair(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_stair(2);
XYcoor_pf(:,1) = XYcoor_pf(:,1) + CoC_stair(1);
XYcoor_pf(:,2) = XYcoor_pf(:,2) + CoC_stair(2);

lengthXYcoor_i = length(XYcoor_i); % 内筒每层节点数
lengthXYcoor_pf = length(XYcoor_pf); % 平台每层节点数
lengthXYcoor_all = lengthXYcoor_i + lengthXYcoor_pf;  % 每层节点数备份

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for j = 1:lengthXYcoor_i % 内部4个点
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
    end
    for j = 1:lengthXYcoor_pf % 外部4个点
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_pf(j,1),XYcoor_pf(j,2),levelZaxis(i));
    end
end
iNO_end = iNO;
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
for i = 1:lengthlevelZaxis	%
    if rem(i,2) ~= 0 % i是奇数层
        iN1 = iNO+1+lengthXYcoor_all*(i-1); % 内筒点1
        iN2 = iN1+3;    % 内筒点4
    else % i是偶数层
        iN1 = iNO+2+lengthXYcoor_all*(i-1); % 内筒点2
        iN2 = iN1+1;    % 内筒点3
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % 梁单元的两个节点号
        ELE_ANGLE, ELE_iSUB);
end

fprintf(fileID,'; 楼梯悬挑梁\n');
iNO = iNO_init; % 初始化iNO
for i = 1:lengthlevelZaxis	%
    if rem(i,2) ~= 0 % i是奇数层
        iN_i(1) = iNO+1+lengthXYcoor_all*(i-1); % 内筒点1
        iN_i(2) = iN_i(1)+3;    % 内筒点4
        iN_pf(1) = iNO+lengthXYcoor_i+1+lengthXYcoor_all*(i-1); % 平台点1
        iN_pf(2) = iN_pf(1)+3;    % 平台点4
    else % i是偶数层
        iN_i(1) = iNO+2+lengthXYcoor_all*(i-1); % 内筒点2
        iN_i(2) = iN_i(1)+1;    % 内筒点3
        iN_pf(1) = iNO+lengthXYcoor_i+2+lengthXYcoor_all*(i-1); % 平台点2
        iN_pf(2) = iN_pf(1)+1;    % 平台点3
    end
    for j = 1:2 % 两根悬挑梁
        iEL = iEL+1;
        iN1 = iN_i(j);
        iN2 = iN_pf(j);
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

fprintf(fileID,'; 楼梯斜梁\n');
iNO = iNO_init; % 初始化iNO
iN_i = zeros(1,lengthXYcoor_i);
iN_table = [1,2; 4,3; 2,1; 3,4]; % 表驱动
for i = 1:(lengthlevelZaxis-1)	% 斜梁，故-1.(逻辑同柱)
    for k = 1:lengthXYcoor_i % 内筒4个点
        iN_i(k) = iNO+k+lengthXYcoor_all*(i-1); % 内筒点1/2/3/4
    end
    for j = 1:2 % 两根斜梁
    if rem(i,2) ~= 0 % i是奇数层
        iN1 = iN_i( iN_table(j,1) ); % 内筒点1/4
        iN2 = iN_i( iN_table(j,2) ) + lengthXYcoor_all;    % 上层内筒点2/3
    else % i是偶数层
        iN1 = iN_i( iN_table(j+2,1) ); % 内筒点2/3
        iN2 = iN_i( iN_table(j+2,2) ) + lengthXYcoor_all;    % 上层内筒点1/4
    end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

% 环形次梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 平台封边梁\n');
ELE_iPRO = 4;
iNO = iNO_init; % 初始化iNO
for i = 1:lengthlevelZaxis	%
    if rem(i,2) ~= 0 % i是奇数层
        iN1 = iNO+lengthXYcoor_i+1+lengthXYcoor_all*(i-1); % 平台点1
        iN2 = iN1+3;    % 平台点4
    else % i是偶数层
        iN1 = iNO+lengthXYcoor_i+2+lengthXYcoor_all*(i-1); % 平台点2
        iN2 = iN1+1;    % 平台点3
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % 梁单元的两个节点号
        ELE_ANGLE, ELE_iSUB);
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
fprintf(fileID,'; 1厚板楼梯斜板\n');
iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis-1) % 由于有斜段，故这里要-1
    if rem(i,2) ~= 0    % 奇数层 1/2+/3+/4
        iN1 = iNO +1 + lengthXYcoor_all*(i-1); % 本层1
        iN2 = iN1 +1 + lengthXYcoor_all; % 上层2
        iN3 = iN2 +1; % 上层3
        iN4 = iN1 +3; % 本层4
    else % 偶数层 1+/2/3/4+
        iN1 = iNO +1 + lengthXYcoor_all*i; % 上层1
        iN2 = iNO +2 + lengthXYcoor_all*(i-1); % 本层2
        iN3 = iN2 +1; % 本层3
        iN4 = iN1 +3; % 上层4
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2, iN3, iN4,...    % 板单元的四个节点号
        ELE_iSUB, ELE_iWID);
end

fprintf(fileID,'; 1厚板楼梯平台板\n');
iNO = iNO_init; % 初始化iNO
for i = 1:lengthlevelZaxis % 平台不需要-1
    if rem(i,2) ~= 0    % 奇数层 1/4/4'/1'
        iN1 = iNO +1 + lengthXYcoor_all*(i-1); % 内筒1
        iN2 = iN1 +3; % 内筒4
        iN3 = iN2 + lengthXYcoor_i; % 平台4
        iN4 = iN1 + lengthXYcoor_i; % 平台1
    else % 偶数层 2/2'/3'/3
        iN1 = iNO +2 + lengthXYcoor_all*(i-1); % 内筒2
        iN2 = iN1 + lengthXYcoor_i; % 平台2
        iN3 = iN2 +1; % 平台3
        iN4 = iN1 +1; % 内筒3
    end
    iEL = iEL+1;
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

fprintf(fileID,'; 1厚板楼梯斜板荷载\n');
iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis-1) % 由于有斜段，故这里要-1
    if rem(i,2) ~= 0    % 奇数层 1/2+/3+/4
        iN1 = iNO +1 + lengthXYcoor_all*(i-1); % 本层1
        iN2 = iN1 +1 + lengthXYcoor_all; % 上层2
        iN3 = iN2 +1; % 上层3
        iN4 = iN1 +3; % 本层4
    else % 偶数层 1+/2/3/4+
        iN1 = iNO +1 + lengthXYcoor_all*i; % 上层1
        iN2 = iNO +2 + lengthXYcoor_all*(i-1); % 本层2
        iN3 = iN2 +1; % 本层3
        iN4 = iN1 +3; % 上层4
    end
    fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
end

fprintf(fileID,'; 1厚板楼梯平台板荷载\n');
iNO = iNO_init; % 初始化iNO
for i = 1:lengthlevelZaxis % 平台不需要-1
    if rem(i,2) ~= 0    % 奇数层 1/4/4'/1'
        iN1 = iNO +1 + lengthXYcoor_all*(i-1); % 内筒1
        iN2 = iN1 +3; % 内筒4
        iN3 = iN2 + lengthXYcoor_i; % 平台4
        iN4 = iN1 + lengthXYcoor_i; % 平台1
    else % 偶数层 2/2'/3'/3
        iN1 = iNO +2 + lengthXYcoor_all*(i-1); % 内筒2
        iN2 = iN1 + lengthXYcoor_i; % 平台2
        iN3 = iN2 +1; % 平台3
        iN4 = iN1 +1; % 内筒3
    end
    fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

end