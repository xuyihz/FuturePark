%% function
% MGT facade
%
% Xu Yi, 28th March 2018
% Xu Yi, 28rd March 2018, revised 全部未改

%%
function [iNO_end, iEL_end] = MGT_facade(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, levelZaxis, levelPstart)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
XYcoor_i = zeros(car_num,2);     % 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o = zeros(car_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。

car_num2pi = 2*pi/car_num;  % speed up

XYcor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % 单个Y型模块内筒一点 X
XYcor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y
XYcor_o_1(1,1) = sqrt(tube_outerR^2 - XYcor_i_1(1,2)^2);        % 单个Y型模块外筒一点 X1 注意外筒16个点并不是等角度等分。
XYcor_o_1(1,2) = XYcor_i_1(1,2);                                % Y1
XYcor_o_1(2,:) = coorMir(XYcor_o_1(1,:), [0,0], XYcor_i_1);     % X2,Y2

for i = 0:(car_num-1)   % 尝试向量化 % 旋转局部角度+整体角度
    [XYcoor_i(i+1,1), XYcoor_i(i+1,2)] = coorTrans(XYcor_i_1(1), XYcor_i_1(2), -car_num2pi*i+Deg_tower);       % 内筒点坐标
    [XYcoor_o(i*2+1,1), XYcoor_o(i*2+1,2)] = coorTrans(XYcor_o_1(1,1), XYcor_o_1(1,2), -car_num2pi*i+Deg_tower); % 外筒点坐标1
    [XYcoor_o(i*2+2,1), XYcoor_o(i*2+2,2)] = coorTrans(XYcor_o_1(2,1), XYcor_o_1(2,2), -car_num2pi*i+Deg_tower); % 外筒点坐标2
end
% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_tower(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_tower(2);
XYcoor_o(:,1) = XYcoor_o(:,1) + CoC_tower(1);
XYcoor_o(:,2) = XYcoor_o(:,2) + CoC_tower(2);

lengthXYcor2 = length(XYcoor_i(:))/2 + length(XYcoor_o(:))/2;  % 每层节点数
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
    for j = 1:car_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+car_num+j+lengthXYcor2*(i-1); % 此行与内筒不同，多了 +car_num
        iN2 = iN1+lengthXYcor2;
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
for i = levelPstart:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1
    for j = 1:car_num	% 每层内筒的节点数
        for k = 1:2 % 一根内筒柱连接两根外筒柱，即两根梁
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcor2*(i-1);
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
        iN1 = iNO+j+lengthXYcor2*(i-1);
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
for i = levelPstart:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1;
    for j = 1:car_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+car_num+j+lengthXYcor2*(i-1); % 此行与内环梁不同，多加了car_num
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
iEL_end = iEL;
fprintf(fileID,'\n');
end