%% function
% MGT tower
%
% Xu Yi, 19th March 2018
% Xu Yi, 23rd April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_tower(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, tube_outerR, levelZaxis, levelPstart1o2, CAR, ~, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
XYcoor_i = zeros(car_num,2);     % 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o = zeros(car_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。

car_num2pi = 2*pi/car_num;  % speed up

XYcoor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % 单个Y型模块内筒一点 X
XYcoor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y
XYcoor_o_1(1,1) = sqrt(tube_outerR^2 - XYcoor_i_1(1,2)^2);        % 单个Y型模块外筒一点 X1 注意外筒16个点并不是等角度等分。
XYcoor_o_1(1,2) = XYcoor_i_1(1,2);                                % Y1
XYcoor_o_1(2,:) = coorMir(XYcoor_o_1(1,:), [0,0], XYcoor_i_1);     % X2,Y2

for i = 0:(car_num-1)   % 尝试向量化 % 旋转局部角度+整体角度
    [XYcoor_i(i+1,:)] = coorTrans(XYcoor_i_1, -car_num2pi*i+Deg_tower);       % 内筒点坐标
    [XYcoor_o(i*2+1,:)] = coorTrans(XYcoor_o_1(1,:), -car_num2pi*i+Deg_tower); % 外筒点坐标1
    [XYcoor_o(i*2+2,:)] = coorTrans(XYcoor_o_1(2,:), -car_num2pi*i+Deg_tower); % 外筒点坐标2
end
% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_tower(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_tower(2);
XYcoor_o(:,1) = XYcoor_o(:,1) + CoC_tower(1);
XYcoor_o(:,2) = XYcoor_o(:,2) + CoC_tower(2);

lengthXYcoor2 = length(XYcoor_i(:))/2 + length(XYcoor_o(:))/2;  % 每层节点数
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for k = 1:2 % 内筒1 外筒2
        for j = 1:car_num
            iNO = iNO+1;
            if k == 1                                           % 节点编号规则：从0度角开始逆时针；先每层内筒，再每层外筒；从下到上。
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
            else
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_o(j*2-1,1),XYcoor_o(j*2-1,2),levelZaxis(i));   % 外筒 X1 & Y1
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_o(j*2,1),XYcoor_o(j*2,2),levelZaxis(i));       % 外筒 X2 & Y2
            end
        end
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
for i = 1:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数
    for j = 1:car_num	% 每层内筒的节点数
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
for i = levelPstart1o2:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数 % levelPstart 第几层开始停车，即下几层开敞
    for j = 1:car_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+car_num+j+lengthXYcoor2*(i-1); % 此行与内筒不同，多了 +car_num
        iN2 = iN1+lengthXYcoor2;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 柱单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) beams
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 横向主梁；iPRO = 3 截面编号3。
fprintf(fileID,'; 横向主梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = levelPstart1o2:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1
    for j = 1:car_num	% 每层内筒的节点数
        for k = 1:2 % 一根内筒柱连接两根外筒柱，即两根梁
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcoor2*(i-1);
            iN2 = iN1-j+car_num+(j-1)*2+k;    % iN1归到内筒第0点后再加car_num后，即为外筒Y型第0点(即内筒最后一点)
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
    for j = 1:car_num	% 每层内筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor2*(i-1);
        if j ~= car_num
            iN2 = iN1+1;
        else % j = car_num 时， 连接的是本环的第一个点，而不是外环的第一个点。
            iN2 = iN1+1-car_num;
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
for i = levelPstart1o2:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1;
    for j = 1:car_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+car_num+j+lengthXYcoor2*(i-1); % 此行与内环梁不同，多加了car_num
        if j ~= car_num*2
            iN2 = iN1+1;
        else % j = car_num*2 时， 连接的是本环的第一个点，而不是上层内环的第一个点。
            iN2 = iN1+1-car_num*2;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) bracings
% iEL = bracings(fileID, iNO_init, iEL, car_num, lengthlevelZaxis, levelPstart, lengthXYcoor2); % 螺旋撑不是必要的设备构件，根据结构计算需要，决定加与不加。

%% ELEMENT(planner) floor
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_floor = iEL;
ELE_TYPE = 'PLATE'; ELE_iMAT = 2; ELE_iSUB = 2; ELE_iWID = 0; % iMAT = 2材料混凝土C30 % iSUB = 2 薄板

% 板厚1；iPRO = 2 截面编号2。
fprintf(fileID,'; 1厚板停车板\n');
ELE_iPRO = 2;
iNO = iNO_init; % 初始化iNO
for i = levelPstart1o2:(lengthlevelZaxis-1) % 此行同外环梁，屋面层先不建即(lengthlevelZaxis-1)
    for j = 1:car_num	% 每层停车数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor2*(i-1);     % 逆时针板四周四个点
        iN2 = iN1-j+1+car_num+(j-1)*2+1;	% iN1归到内筒第一点后再加car_num后，即为外筒第一点(即Y型第一点，实际起点应再+1，即为Y型第二点)
        if j ~= car_num
            iN3 = iN2+1;
            iN4 = iN1+1;
        else % j = car_num 时， 连接的是这层的第一个点，而不是上层的第一个点。
            iN3 = iN2+1-car_num*2;
            iN4 = iN1+1-car_num;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2, iN3, iN4,...    % 板单元的四个节点号
            ELE_iSUB, ELE_iWID);
    end
end
iEL_end = iEL;
fprintf(fileID,'\n');

%% FLOORLOAD
fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');

LTNAME = CAR; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0;
DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';

iNO = iNO_init; % 初始化iNO
for i = levelPstart1o2:(lengthlevelZaxis-1) % 此行同1厚板停车板，即同外环梁，屋面层先不建即(lengthlevelZaxis-1)
    for j = 1:car_num	% 每层停车数
        iN1 = iNO+j+lengthXYcoor2*(i-1); % 逆时针板四周四个点
        iN2 = iN1-j+1+car_num+(j-1)*2+1;	% iN1归到内筒第一点后再加car_num后，即为外筒第一点(即Y型第一点，实际起点应再+1，即为Y型第二点)
        if j ~= car_num
            iN3 = iN2+1;
            iN4 = iN1+1;
        else % j = car_num 时， 连接的是这层的第一个点，而不是上层的第一个点。
            iN3 = iN2+1-car_num*2;
            iN4 = iN1+1-car_num;
        end
        fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
    end
end
fprintf(fileID,'\n');

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % 初始化iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+car_num);
CONSTRAINT = 111111; % 6个自由度全约束
fprintf(fileID,'   %s, %d, \n',...
                NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end

%% ELEMENT(frame) bracings "8->16"后未修改
% function  iEL = bracings(fileID, iNO_init, iEL, car_num, lengthlevelZaxis, levelPstart, lengthXYcoor2)
% fprintf(fileID,'*ELEMENT    ; Elements\n');
% fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');
% 
% % iEL_init_bracing = iEL;
% ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345
% 
% % 斜撑；iPRO = 5 截面编号5。
% fprintf(fileID,'; 螺旋撑\n');
% ELE_iPRO = 5;
% iNO = iNO_init; % 初始化iNO
% for i = levelPstart:2:(lengthlevelZaxis-2)	% 由于螺旋撑为每两层一根，故为间隔2； 此处与柱梁都不同，因两层一撑，故要-2
%     for j = 1:car_num	% 每层外筒的节点数
%         iEL = iEL+1;
%         iN1 = iNO+(j+car_num)+lengthXYcoor2*(i-1); % 此行与柱单元相同
%         if j ~= car_num
%             iN2 = iN1+1+lengthXYcoor2*2;
%         else % j = car_num 时， 连接的是这层外环的第一个点，而不是上层内环的第一个点。
%             iN2 = iN1+1+lengthXYcoor2*2-car_num;
%         end
%         fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
%             iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
%             iN1, iN2,...    % 梁单元的两个节点号
%             ELE_ANGLE, ELE_iSUB);
%     end
% end
% fprintf(fileID,'\n');
% end