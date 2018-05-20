%% function
% MGT towerC ramp Misc.
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_Misc(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, ~, levelZaxis, levelPstart, Roof_boundary, ~)    %注意这里的levelPstart是1x2数组
%% NODE 定义商业层和屋面层上坡道投影处的节点
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

% 由于坡道每层外挑都有变化，故需要循环(或向量化)。(暂时未考虑外挑变化) % 从17600开始往下、往上建。
ramp_inner_R = 12200; % 坡道内圈半径
ramp_outer_R = 14400; % 坡道外圈半径
aisle_width = 1800; % 走道宽度

iNO_init = iNO;
% levelPstart1 = levelPstart(1);
levelPstart2 = levelPstart(2);
% XYcoor_i = zeros(car_num,2);	% 内筒XoY坐标第1(X)、2(Y)列。
% XYcoor_o = zeros(car_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。
ramp_point = 1;
XYcoor_ramp_i = zeros(ramp_point,car_num*2,2);	% 坡道内侧XoY坐标第1(X)、2(Y)列。
XYcoor_ramp_o = zeros(ramp_point,car_num*2,2);	% 坡道外侧XoY坐标第1(X)、2(Y)列。
XYcoor_aisle = zeros(ramp_point,car_num*2,2);	% 走道外侧XoY坐标第1(X)、2(Y)列。(走道内侧即为坡道外侧)
circle_num = 3; % 以上三个圈

car_num2pi = 2*pi/car_num;  % speed up
Deg_tower = Deg_tower + pi/2; % 调整坡道定位

XYcoor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % 单个Y型模块内筒一点 X % Y型模块看Parts文件夹内单筒文件
XYcoor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y

XYcoor_ramp_i_1(1,1) = sqrt(ramp_inner_R^2 - XYcoor_i_1(1,2)^2);        % Y型坡道内侧上一点 X1 注意16个点并不是等角度等分。
XYcoor_ramp_i_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_ramp_i_1(2,:) = coorMir(XYcoor_ramp_i_1(1,:), [0,0], XYcoor_i_1);% 对Y的主干线镜像，得到Y型模块下面分支的X2,Y2
XYcoor_ramp_o_1(1,1) = sqrt(ramp_outer_R^2 - XYcoor_i_1(1,2)^2);        % Y型坡道外侧上一点 X1 注意16个点并不是等角度等分。
XYcoor_ramp_o_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_ramp_o_1(2,:) = coorMir(XYcoor_ramp_o_1(1,:), [0,0], XYcoor_i_1);% 对Y的主干线镜像，得到Y型模块下面分支的X2,Y2
XYcoor_aisle_1(1,1) = sqrt((ramp_outer_R+aisle_width)^2 - XYcoor_i_1(1,2)^2);        % Y型走道外侧上一点 X1 注意16个点并不是等角度等分。
XYcoor_aisle_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_aisle_1(2,:) = coorMir(XYcoor_aisle_1(1,:), [0,0], XYcoor_i_1);% 对Y的主干线镜像，得到Y型模块下面分支的X2,Y2
for j = 1:ramp_point % 暂定1，由于坡道悬挑长度和高度有关，后期还会更改。
    for i = 0:(car_num-1)   % 尝试向量化 % 旋转局部角度+整体角度
        [XYcoor_ramp_i(j,i*2+1,:)] = coorTrans(XYcoor_ramp_i_1(1,:), -car_num2pi*i+Deg_tower); % 坡道内侧点坐标1
        [XYcoor_ramp_i(j,i*2+2,:)] = coorTrans(XYcoor_ramp_i_1(2,:), -car_num2pi*i+Deg_tower); % 坡道内侧点坐标2
        [XYcoor_ramp_o(j,i*2+1,:)] = coorTrans(XYcoor_ramp_o_1(1,:), -car_num2pi*i+Deg_tower); % 坡道外侧点坐标1
        [XYcoor_ramp_o(j,i*2+2,:)] = coorTrans(XYcoor_ramp_o_1(2,:), -car_num2pi*i+Deg_tower); % 坡道外侧点坐标2
        [XYcoor_aisle(j,i*2+1,:)] = coorTrans(XYcoor_aisle_1(1,:), -car_num2pi*i+Deg_tower); % 坡道外侧点坐标1
        [XYcoor_aisle(j,i*2+2,:)] = coorTrans(XYcoor_aisle_1(2,:), -car_num2pi*i+Deg_tower); % 坡道外侧点坐标2
    end
end
% 局部坐标系 转换至 整体坐标系
XYcoor_ramp_i(:,:,1) = XYcoor_ramp_i(:,:,1) + CoC_tower(1);
XYcoor_ramp_i(:,:,2) = XYcoor_ramp_i(:,:,2) + CoC_tower(2);
XYcoor_ramp_o(:,:,1) = XYcoor_ramp_o(:,:,1) + CoC_tower(1);
XYcoor_ramp_o(:,:,2) = XYcoor_ramp_o(:,:,2) + CoC_tower(2);
XYcoor_aisle(:,:,1) = XYcoor_aisle(:,:,1) + CoC_tower(1);
XYcoor_aisle(:,:,2) = XYcoor_aisle(:,:,2) + CoC_tower(2);

[~,lengthXYcoor_ramp,~] = size(XYcoor_ramp_i);  % 坡道每层节点数
for k = 1:2
    if k == 1
        levelZ = levelZaxis(levelPstart2);
    else
        levelZ = levelZaxis(end);
    end
    i = ramp_point;
    for j = 1:lengthXYcoor_ramp  % 商业层、屋面层 坡道x悬臂梁投影节点
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_ramp_i(i,j,1),XYcoor_ramp_i(i,j,2),levelZ);   % 外筒 X & Y
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_ramp_o(i,j,1),XYcoor_ramp_o(i,j,2),levelZ);   % 外筒 X & Y
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_aisle(i,j,1),XYcoor_aisle(i,j,2),levelZ);   % 外筒 X & Y
    end
end
fprintf(fileID,'\n');

iNO_init2 = iNO;
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');
% 定义屋面边线
% Roof_boundary = [30833,4750; 5248,41010; 45186,54754; 77970,61552; 87304,37382; 117443,7243; 114951,4750]; % 外边线定位点，从左下角点起，顺时针定位点
levelZ = levelZaxis(end);
for i = 1:length(Roof_boundary)
    iNO = iNO+1;
    fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
        iNO,Roof_boundary(i,1),Roof_boundary(i,2),levelZ);   % 外筒 X & Y
end

iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns 

%% ELEMENT(frame) beams 定义商业层和屋面层上坡道投影处的环梁
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 环形次梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 环形次梁\n');
ELE_iPRO = 4;
iNO = iNO_init; % 初始化iNO
% 外环梁
fprintf(fileID,';   坡道环梁\n');
for k = 1:2 % 商业层、屋面层
    for i = 1:lengthXYcoor_ramp % 坡道投影一圈的节点数
        for j = 1:circle_num % 坡道内/外环梁/走道外环梁
            iEL = iEL+1;
            iN1 = iNO+j+(i-1)*circle_num+(k-1)*lengthXYcoor_ramp*circle_num;	% 外筒上的点/坡道内上的点
            if i == lengthXYcoor_ramp
                iN2 = iN1+circle_num-lengthXYcoor_ramp*circle_num;
            else
                iN2 = iN1+circle_num;    % 坡道内上的点/坡道外/走道外上的点
            end
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 梁单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        end
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');
% 屋面环形边梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 屋面环形边梁\n');
ELE_iPRO = 4;
iNO = iNO_init2; % 初始化iNO
for i = 1:length(Roof_boundary)
    iEL = iEL+1;
    iN1 = iNO+i;	% 一字长蛇阵
    if i == length(Roof_boundary)
        iN2 = iNO+1;
    else
        iN2 = iN1+1;    % 
    end
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % 梁单元的两个节点号
        ELE_ANGLE, ELE_iSUB);
end
fprintf(fileID,'\n');

iEL_end = iEL;
end
