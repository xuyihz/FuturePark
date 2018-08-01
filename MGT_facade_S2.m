%% function
% MGT tower facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_facade_S2(fileID, iNO, iEL, column_num, CoC_tower, Deg_tower, facade_tower_R, tube_innerR, levelZaxis, levelPstart, iNO_towerS_init)    %注意这里的levelPstart是1x2数组
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
levelPstart1 = levelPstart(1); levelPstart2 = levelPstart(2);
lengthlevelZaxis = length(levelZaxis(:));
XYcoor_i = zeros(column_num,2);    % 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o3 = zeros(lengthlevelZaxis,column_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。注意这里是三维数组。(Z方向幕墙有变化)

car_num2pi = 2*pi/column_num;  % speed up

XYcoor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % 单个Y型模块内筒一点 X
XYcoor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y

for j = levelPstart1:lengthlevelZaxis
    XYcoor_o_1(1,1) = sqrt(facade_tower_R(j)^2 - XYcoor_i_1(1,2)^2);        % Y型幕墙上一点 X1 注意幕墙16个点并不是等角度等分。
    XYcoor_o_1(1,2) = XYcoor_i_1(1,2);                                % Y1
    XYcoor_o_1(2,:) = coorMir(XYcoor_o_1(1,:), [0,0], XYcoor_i_1);     % X2,Y2
    for i = 0:(column_num-1)   % 尝试向量化 % 旋转局部角度+整体角度
        [XYcoor_i(i+1,:)] = coorTrans(XYcoor_i_1, -car_num2pi*i+Deg_tower);       % 内筒点坐标 % 已旋转整体角度
        [XYcoor_o3(j,i*2+1,:)] = coorTrans(XYcoor_o_1(1,:), -car_num2pi*i+Deg_tower); % 外筒点坐标1
        [XYcoor_o3(j,i*2+2,:)] = coorTrans(XYcoor_o_1(2,:), -car_num2pi*i+Deg_tower); % 外筒点坐标2
    end
end
% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_tower(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_tower(2);
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_tower(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_tower(2);

lengthXYcoor2 = length(XYcoor_i(:))/2 + column_num*2;  % 每层节点数备份
lengthXYcoor_f = column_num*2;  % 幕墙每层节点数

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for j = 1:(column_num*2)
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));   % 外筒 X & Y
    end
end
iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_colu = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 外筒柱；iPRO = 1 截面编号1。
fprintf(fileID,'; 幕墙柱(类似外筒柱)\n');
ELE_iPRO = 1;
iNO = iNO_init; % 初始化iNO
for i = levelPstart1:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数 % levelPstart 第几层开始停车，即下几层开敞
    for j = 1:column_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_f*(i-1); % 此行幕墙与另外塔的不同，因幕墙节点为单独编号，故一层只有lengthXYcoor_f个节点
        iN2 = iN1+lengthXYcoor_f;
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
iNO_towerS = iNO_towerS_init; % 初始化iNO_towerS
for i = levelPstart1:(levelPstart2-1)	% 不停车层，梁从内筒伸出
    for j = 1:column_num	% 每层内筒的节点数
        for k = 1:2 % 一根内筒柱连接两根外筒柱，即两根梁
            iEL = iEL+1;
            iN1 = iNO_towerS+j+lengthXYcoor2*(i-1); % 此行为定位梁在塔楼的节点(内筒)
            iN2 = iNO+lengthXYcoor_f*(i-1)+(j-1)*2+k;    % 归到幕墙外筒第0点后，再定位到具体点
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 梁单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        end
    end
end
for i = levelPstart2:lengthlevelZaxis	% 停车层，梁从外筒伸出
    for j = 1:column_num	% 每层内筒的节点数
        for k = 1:2 % 两根外筒柱，即两根梁
            iEL = iEL+1;
            iN1 = iNO_towerS+column_num+lengthXYcoor2*(i-1)+(j-1)*2+k; % 此行为定位梁在塔楼的节点(外筒)
            iN2 = iNO+lengthXYcoor_f*(i-1)+(j-1)*2+k;    % 归到幕墙外筒第0点后，再定位到具体点
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
for i = levelPstart1:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1;
    for j = 1:column_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_f*(i-1); % 
        if j ~= column_num*2
            iN2 = iN1+1;
        else % j = car_num*2 时， 连接的是本环的第一个点，而不是上层内环的第一个点。
            iN2 = iN1+1-column_num*2;
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