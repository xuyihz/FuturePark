%% function
% MGT elevator
%
% Xu Yi, 24th April 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_elevator_platform(fileID, iNO, iEL, CoC_elevator, Deg_elevator, levelZaxis, elevatorColu_num, elevatorR, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis = length(levelZaxis(:));

XYcoor_i = zeros(elevatorColu_num,2);   % 8个点(含电梯中间不落柱的1个点和圆心点)内筒XoY坐标第1(X)、2(Y)列。

% 内筒点
elevatorXYtemp = elevatorR/sqrt(2);
elevatorXY = [elevatorXYtemp,elevatorXYtemp; -elevatorXYtemp,elevatorXYtemp;...
    -elevatorXYtemp,-elevatorXYtemp; elevatorXYtemp,-elevatorXYtemp;...
    elevatorXYtemp,350; -elevatorXYtemp,350; 0,elevatorXYtemp; 0,350]; % 原始坐标，未旋转，未转到整体坐标系
for i = 1:elevatorColu_num   % 尝试向量化
    XYcoor_i(i,:) = coorTrans(elevatorXY(i,:), Deg_elevator); % 内筒
end

% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_elevator(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_elevator(2);

lengthXYcoor_i = length(XYcoor_i); % 内筒每层节点数
lengthXYcoor_all = lengthXYcoor_i;  % 每层节点数备份

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    % 内筒 从下到上。
    for j = 1:lengthXYcoor_i % 内部7个柱子(有一个不落柱)
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
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
fprintf(fileID,'; 电梯筒楼梯斜梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
iN_i = zeros(1,lengthXYcoor_i);
iN_table = [5,6; 4,3; 6,5; 3,4]; % 表驱动 % 每行为梁两端节点号
for i = 1:(lengthlevelZaxis-1)	% 斜梁，故-1.(逻辑同柱)
    for k = 1:lengthXYcoor_i % 内筒8个点
        iN_i(k) = iNO+k+lengthXYcoor_all*(i-1); % 内筒点1~8
    end
    for j = 1:2 % 两根斜梁
        if rem(i,2) ~= 0 % i是奇数层
            iN1 = iN_i( iN_table(j,1) ); % 内筒点5/4
            iN2 = iN_i( iN_table(j,2) ) + lengthXYcoor_all;    % 上层内筒点6/3
        else % i是偶数层
            iN1 = iN_i( iN_table(j+2,1) ); % 内筒点6/3
            iN2 = iN_i( iN_table(j+2,2) ) + lengthXYcoor_all;    % 上层内筒点5/4
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

fprintf(fileID,'; 电梯筒楼梯封边梁\n');
iNO = iNO_init; % 初始化iNO
for i = 1:lengthlevelZaxis	%
    if rem(i,2) ~= 0 % i是奇数层
        iN1 = iNO+4+lengthXYcoor_all*(i-1); % 内筒点4
        iN2 = iN1+1;    % 内筒点5
    else % i是偶数层
        iN1 = iNO+3+lengthXYcoor_all*(i-1); % 内筒点3
        iN2 = iN1+3;    % 内筒点6
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
    if rem(i,2) ~= 0    % 奇数层 5/6+/3+/4
        iN1 = iNO +5 + lengthXYcoor_all*(i-1); % 本层5
        iN2 = iN1 +1 + lengthXYcoor_all; % 上层6
        iN3 = iN2 -3; % 上层3
        iN4 = iN1 -1; % 本层4
    else % 偶数层 5+/6/3/4+
        iN1 = iNO +5 + lengthXYcoor_all*i; % 上层5
        iN2 = iNO +6 + lengthXYcoor_all*(i-1); % 本层6
        iN3 = iN2 -3; % 本层3
        iN4 = iN1 -1; % 上层4
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2, iN3, iN4,...    % 板单元的四个节点号
        ELE_iSUB, ELE_iWID);
end
fprintf(fileID,'\n');

%% FLOORLOAD 楼面荷载加不上
fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');

LTNAME = ROOF; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0; % 楼梯荷载 目前采用了屋面荷载4/3.5，板厚0.因每2.1一块板，故可能和4.2一层7/3.5差不多。(待复核)
DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';

fprintf(fileID,'; 1厚板楼梯斜板荷载\n');
iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis-1) % 由于有斜段，故这里要-1
    if rem(i,2) ~= 0    % 奇数层 5/6+/3+/4
        iN1 = iNO +5 + lengthXYcoor_all*(i-1); % 本层5
        iN2 = iN1 +1 + lengthXYcoor_all; % 上层6
        iN3 = iN2 -3; % 上层3
        iN4 = iN1 -1; % 本层4
    else % 偶数层 5+/6/3/4+
        iN1 = iNO +5 + lengthXYcoor_all*i; % 上层5
        iN2 = iNO +6 + lengthXYcoor_all*(i-1); % 本层6
        iN3 = iN2 -3; % 本层3
        iN4 = iN1 -1; % 上层4
    end
    fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

end