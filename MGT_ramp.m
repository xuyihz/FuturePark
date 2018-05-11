%% function
% MGT towerC ramp
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_ramp(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, tube_outerR, levelZaxis, levelPstart, iNO_towerS_init)    %注意这里的levelPstart是1x2数组
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

% 由于坡道每层外挑都有变化，故需要循环(或向量化)。(暂时未考虑外挑变化) % 从17600开始往下、往上建。
ramp_inner_R = 12200; % 坡道内圈半径
ramp_outer_R = 14400; % 坡道外圈半径

iNO_init = iNO;
levelPstart1 = levelPstart(1); levelPstart2 = levelPstart(2);
lengthlevelZaxis = length(levelZaxis(:));
% XYcoor_i = zeros(car_num,2);	% 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o = zeros(car_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。
XYcoor_ramp_i = zeros(1,car_num*2,2);	% 坡道内侧XoY坐标第1(X)、2(Y)列。
XYcoor_ramp_o = zeros(1,car_num*2,2);	% 坡道外侧XoY坐标第1(X)、2(Y)列。

car_num2pi = 2*pi/car_num;  % speed up

XYcoor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % 单个Y型模块内筒一点 X % Y型模块看Parts文件夹内单筒文件
XYcoor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y

XYcoor_o_1(1,1) = sqrt(tube_outerR^2 - XYcoor_i_1(1,2)^2);        % 单个Y型模块外筒一点 X1 注意外筒16个点并不是等角度等分。
XYcoor_o_1(1,2) = XYcoor_i_1(1,2);                                % Y1
XYcoor_o_1(2,:) = coorMir(XYcoor_o_1(1,:), [0,0], XYcoor_i_1);     % X2,Y2
for i = 0:(car_num-1)   % 尝试向量化 % 旋转局部角度+整体角度
    [XYcoor_o(i*2+1,1), XYcoor_o(i*2+1,2)] = coorTrans(XYcoor_o_1(1,1), XYcoor_o_1(1,2), -car_num2pi*i+Deg_tower); % 外筒点坐标1
    [XYcoor_o(i*2+2,1), XYcoor_o(i*2+2,2)] = coorTrans(XYcoor_o_1(2,1), XYcoor_o_1(2,2), -car_num2pi*i+Deg_tower); % 外筒点坐标2
end

XYcoor_ramp_i_1(1,1) = sqrt(ramp_inner_R^2 - XYcoor_i_1(1,2)^2);        % Y型坡道内侧上一点 X1 注意16个点并不是等角度等分。
XYcoor_ramp_i_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_ramp_i_1(2,:) = coorMir(XYcoor_ramp_i_1(1,:), [0,0], XYcoor_i_1);% 对Y的主干线镜像，得到Y型模块下面分支的X2,Y2
XYcoor_ramp_o_1(1,1) = sqrt(ramp_outer_R^2 - XYcoor_i_1(1,2)^2);        % Y型坡道外侧上一点 X1 注意16个点并不是等角度等分。
XYcoor_ramp_o_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_ramp_o_1(2,:) = coorMir(XYcoor_ramp_o_1(1,:), [0,0], XYcoor_i_1);% 对Y的主干线镜像，得到Y型模块下面分支的X2,Y2
for j = 1 % 暂定1，由于坡道悬挑长度和高度有关，后期还会更改。
    for i = 0:(car_num-1)   % 尝试向量化 % 旋转局部角度+整体角度
        [XYcoor_ramp_i(j,i*2+1,1), XYcoor_ramp_i(j,i*2+1,2)] = coorTrans(XYcoor_ramp_i_1(1,1), XYcoor_ramp_i_1(1,2), -car_num2pi*i+Deg_tower); % 坡道内侧点坐标1
        [XYcoor_ramp_i(j,i*2+2,1), XYcoor_ramp_i(j,i*2+2,2)] = coorTrans(XYcoor_ramp_i_1(2,1), XYcoor_ramp_i_1(2,2), -car_num2pi*i+Deg_tower); % 坡道内侧点坐标2
        [XYcoor_ramp_o(j,i*2+1,1), XYcoor_ramp_o(j,i*2+1,2)] = coorTrans(XYcoor_ramp_o_1(1,1), XYcoor_ramp_o_1(1,2), -car_num2pi*i+Deg_tower); % 坡道外侧点坐标1
        [XYcoor_ramp_o(j,i*2+2,1), XYcoor_ramp_o(j,i*2+2,2)] = coorTrans(XYcoor_ramp_o_1(2,1), XYcoor_ramp_o_1(2,2), -car_num2pi*i+Deg_tower); % 坡道外侧点坐标2
    end
end
% 局部坐标系 转换至 整体坐标系
% XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_tower(1);
% XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_tower(2);
XYcoor_o(:,1) = XYcoor_o(:,1) + CoC_tower(1);
XYcoor_o(:,2) = XYcoor_o(:,2) + CoC_tower(2);
XYcoor_ramp_i(:,:,1) = XYcoor_ramp_i(:,:,1) + CoC_tower(1);
XYcoor_ramp_i(:,:,2) = XYcoor_ramp_i(:,:,2) + CoC_tower(2);
XYcoor_ramp_o(:,:,1) = XYcoor_ramp_o(:,:,1) + CoC_tower(1);
XYcoor_ramp_o(:,:,2) = XYcoor_ramp_o(:,:,2) + CoC_tower(2);

[~,lengthXYcoor_ramp,~] = size(XYcoor_ramp_i);  % 坡道每层节点数
% 坡道悬挑梁标高的确定，暂定一圈坡道6000高
!theta_temp = arctan ( (XYcoor_ramp_i_1(1,2)+XYcoor_ramp_o_1(1,2)) / (XYcoor_ramp_i_1(1,1)+XYcoor_ramp_o_1(1,1)) );


for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for j = 1:lengthXYcoor_ramp
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_o(j,1),XYcoor_o(j,2),levelZaxis(i));   % 外筒 X & Y
    end
    for j = 1:lengthXYcoor_ramp
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_ramp_i(i,j,1),XYcoor_ramp_i(i,j,2),levelZaxis(i));   % 外筒 X & Y
    end
    for j = 1:lengthXYcoor_ramp
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_ramp_o(i,j,1),XYcoor_ramp_o(i,j,2),levelZaxis(i));   % 外筒 X & Y
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
    for j = 1:car_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_ramp*(i-1); % 此行幕墙与另外塔的不同，因幕墙节点为单独编号，故一层只有lengthXYcoor_f个节点
        iN2 = iN1+lengthXYcoor_ramp;
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
    for j = 1:car_num	% 每层内筒的节点数
        for k = 1:2 % 一根内筒柱连接两根外筒柱，即两根梁
            iEL = iEL+1;
            iN1 = iNO_towerS+j+lengthXYcoor2*(i-1); % 此行为定位梁在塔楼的节点(内筒)
            iN2 = iNO+lengthXYcoor_ramp*(i-1)+(j-1)*2+k;    % 归到幕墙外筒第0点后，再定位到具体点
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 梁单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        end
    end
end
for i = levelPstart2:lengthlevelZaxis	% 停车层，梁从外筒伸出
    for j = 1:car_num	% 每层内筒的节点数
        for k = 1:2 % 两根外筒柱，即两根梁
            iEL = iEL+1;
            iN1 = iNO_towerS+car_num+lengthXYcoor2*(i-1)+(j-1)*2+k; % 此行为定位梁在塔楼的节点(外筒)
            iN2 = iNO+lengthXYcoor_ramp*(i-1)+(j-1)*2+k;    % 归到幕墙外筒第0点后，再定位到具体点
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
    for j = 1:car_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_ramp*(i-1); % 
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