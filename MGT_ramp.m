%% function
% MGT towerC ramp
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_ramp(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, tube_outerR, levelZaxis, levelPstart, ~)    %注意这里的levelPstart是1x2数组
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

% 由于坡道每层外挑都有变化，故需要循环(或向量化)。(暂时未考虑外挑变化) % 从17600开始往下、往上建。
ramp_inner_R = 12200; % 坡道内圈半径
ramp_outer_R = 14400; % 坡道外圈半径

iNO_init = iNO;
levelPstart1 = levelPstart(1); levelPstart2 = levelPstart(2);
% XYcoor_i = zeros(car_num,2);	% 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o = zeros(car_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。
ramp_point = 1;
XYcoor_ramp_i = zeros(ramp_point,car_num*2,2);	% 坡道内侧XoY坐标第1(X)、2(Y)列。
XYcoor_ramp_o = zeros(ramp_point,car_num*2,2);	% 坡道外侧XoY坐标第1(X)、2(Y)列。

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
for j = 1:ramp_point % 暂定1，由于坡道悬挑长度和高度有关，后期还会更改。
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
% 坡道悬挑梁Z向标高的确定，暂定一圈坡道6000高
ramp_height_1 = 6000;
theta_temp = atan ( (XYcoor_ramp_i_1(1,2)+XYcoor_ramp_o_1(1,2)) / (XYcoor_ramp_i_1(1,1)+XYcoor_ramp_o_1(1,1)) ); %取悬臂梁在坡道中心线交点的角度
theta_para = theta_temp*2; % 平行悬臂梁parallel与坡道中心线交点的夹角
theta_Y = 2*pi/8 - theta_para; % Y型悬臂梁与坡道中心线交点的夹角
ramp_Z_sample = zeros(1,lengthXYcoor_ramp);
for i = 2:lengthXYcoor_ramp % i-1时，ramp相对标高为0
    if rem(i,2) == 0 % i是偶数
        ramp_beam_deg = theta_Y;
    else
        ramp_beam_deg = theta_para;
    end
    ramp_Z_sample = ramp_Z_sample + [zeros(1,i-1),ones(1,lengthXYcoor_ramp-i+1)*ramp_beam_deg/(2*pi)*ramp_height_1];
end
rampParkArea = levelZaxis(end-2) - levelZaxis(levelPstart1);
length_rampParkArea = fix(rampParkArea/ramp_height_1)*lengthXYcoor_ramp; % fix向零取整
ramp_compare = rem(rampParkArea,ramp_height_1) >= ramp_Z_sample;
for i=1:lengthXYcoor_ramp
    if ramp_compare(i) == 0
        length_rampParkArea = length_rampParkArea + i-1; % 旋转坡道商业层(含)以下至最下停车层的节点数。
        break
    end
end

for jj = 1:length_rampParkArea  % 从商业层往下由内向外顺时针定义
    levelZ = levelZaxis(end-2)-ramp_height_1*fix(jj/lengthXYcoor_ramp)-ramp_Z_sample(rem(jj,lengthXYcoor_ramp)+1);
    i = ramp_point;
    j = lengthXYcoor_ramp - rem(jj,lengthXYcoor_ramp);
    iNO = iNO+1;
    fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
        iNO,XYcoor_o(j,1),XYcoor_o(j,2),levelZ);   % 外筒 X & Y
    iNO = iNO+1;
    fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
        iNO,XYcoor_ramp_i(i,j,1),XYcoor_ramp_i(i,j,2),levelZ);   % 外筒 X & Y
    iNO = iNO+1;
    fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
        iNO,XYcoor_ramp_o(i,j,1),XYcoor_ramp_o(i,j,2),levelZ);   % 外筒 X & Y
end
iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns 悬臂梁没有柱
% fprintf(fileID,'*ELEMENT    ; Elements\n');
% fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');
% 
% % iEL_init_colu = iEL;
% ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345
% 
% % 外筒柱；iPRO = 1 截面编号1。
% fprintf(fileID,'; 幕墙柱(类似外筒柱)\n');
% ELE_iPRO = 1;
% iNO = iNO_init; % 初始化iNO
% for i = levelPstart1:(length_rampParkArea-1)	% length(A(:)) A向量元素个数 % levelPstart 第几层开始停车，即下几层开敞
%     for j = 1:car_num*2	% 每层外筒的节点数
%         iEL = iEL+1;
%         iN1 = iNO+j+lengthXYcoor_ramp*(i-1); % 此行幕墙与另外塔的不同，因幕墙节点为单独编号，故一层只有lengthXYcoor_f个节点
%         iN2 = iN1+lengthXYcoor_ramp;
%         fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
%             iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
%             iN1, iN2,...    % 柱单元的两个节点号
%             ELE_ANGLE, ELE_iSUB);
%     end
% end
% fprintf(fileID,'\n');

%% ELEMENT(frame) beams
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 横向主梁；iPRO = 3 截面编号3。
fprintf(fileID,'; 悬臂梁\n');
ELE_iPRO = 3;
iNO = iNO_init; % 初始化iNO
for i = 1:length_rampParkArea	% 坡道悬臂梁广义层数
    for j = 1:2 % 悬臂梁分两段
        iEL = iEL+1;
        iN1 = iNO+j+(i-1)*3;	% 外筒上的点/坡道内上的点
        iN2 = iN1+1;    % 坡道内上的点/坡道外上的点
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end

% 环形次梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 环形次梁\n');
ELE_iPRO = 4;
iNO = iNO_init; % 初始化iNO
% 外环梁
fprintf(fileID,';   坡道环梁\n');
for i = 1:(length_rampParkArea-1)	% 坡道悬臂梁广义层数
    for j = 1:2 % 坡道内/外环梁
        iEL = iEL+1;
        iN1 = iNO+(j+1)+(i-1)*3;	% 外筒上的点/坡道内上的点
        iN2 = iN1+3;    % 坡道内上的点/坡道外上的点
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
iEL_end = iEL;
fprintf(fileID,'\n');
end