%% function
% MGT stair facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_facade_stair(fileID, iNO, iEL, stairColu_num, CoC_stair, Deg_stair, stairL, stairW, levelZaxis, levelPstart1, iNO_stair_init, tower_num)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

% 由于幕墙每层外挑都有变化，故需要循环(或向量化)。与层数有关，与控制点定位有关，先暂定固定的数值(目前共19层(含半平台)，从7层起到19层)，后需参数化。
switch tower_num % 塔标号
    case {5, 6}
        facade_stair_R = [zeros(6,1); 6820; 5839; 5162; 4734; 4529; 4534; 4750; 5189; 5880; 6878; 8473; 11040; 16200];
    case {7, 9}
        facade_stair_R = [zeros(6,1); 7594; 6144; 5151; 4495; 4120; 4000; 4128; 4511; 5176; 6182; 7830; 10527; 16000];
    case 8
        facade_stair_R = [zeros(6,1); 9769; 8319; 7326; 6669; 6295; 6175; 6302; 6685; 7351; 8357; 10005; 12702; 18175];
    case 10
        facade_stair_R = [zeros(6,1); 6927; 5477; 4483; 3827; 3452; 3333; 3460; 3843; 4509; 5514; 7163; 9860; 15333];
end

iNO_init = iNO;
lengthlevelZaxis = length(levelZaxis(:));
XYcoor_i = zeros(stairColu_num,2);    % 内筒XoY坐标第1(X)、2(Y)列。
XYcoor_o3 = zeros(lengthlevelZaxis,stairColu_num*2,2);	% 外筒XoY坐标第1(X)、2(Y)列。注意这里是三维数组。(Z方向幕墙有变化)

stairXY = [stairL/2, stairW/2; -stairL/2, stairW/2; -stairL/2, -stairW/2; stairL/2, -stairW/2]; % 原始坐标，未旋转，未转到整体坐标系
for i = 1:stairColu_num   % 尝试向量化
    [XYcoor_i(i,1), XYcoor_i(i,2)] = coorTrans(stairXY(i,1), stairXY(i,2), Deg_stair); % 内筒
end
% 幕墙外筒
for j = levelPstart1:lengthlevelZaxis
    XYcoor_o_x = sqrt(facade_stair_R(j)^2 - (stairW/2)^2);  % 幕墙外圈点当y=±stairW/2时，x的绝对值
    XYcoor_o_y = sqrt(facade_stair_R(j)^2 - (stairL/2)^2);  % 幕墙外圈点当x=±stairL/2时，y的绝对值
    stairXY2 = [XYcoor_o_x, stairW/2; stairL/2, XYcoor_o_y; -stairL/2, XYcoor_o_y; -XYcoor_o_x, stairW/2;...
        -XYcoor_o_x, -stairW/2; -stairL/2, -XYcoor_o_y; stairL/2, -XYcoor_o_y; XYcoor_o_x, -stairW/2];  % 幕墙外筒各点坐标
    for i = 1:stairColu_num*2   % 尝试向量化
        [XYcoor_o3(j,i,1), XYcoor_o3(j,i,2)] = coorTrans(stairXY2(i,1), stairXY2(i,2), Deg_stair); % 幕墙外筒点坐标 % 已旋转整体角度
    end
end

% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_stair(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_stair(2);
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_stair(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_stair(2);

lengthXYcoor2 = length(XYcoor_i(:))/2 + stairColu_num*2;  % 每层节点数备份
lengthXYcoor_f = stairColu_num*2;  % 幕墙每层节点数

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for j = 1:(stairColu_num*2)
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));   % 外筒 X1 & Y1
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
    for j = 1:stairColu_num*2	% 每层外筒的节点数
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
iNO_stair = iNO_stair_init; % 初始化iNO_stair
for i = levelPstart1:lengthlevelZaxis	% 楼梯都是从外筒伸出
    for j = 1:stairColu_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO_stair+(stairColu_num+j)+lengthXYcoor2*(i-1); % 此行为定位梁在塔楼的节点(外筒)
        iN2 = iNO+lengthXYcoor_f*(i-1)+j;    % 归到幕墙外筒第0点后，再定位到具体点
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
fprintf(fileID,';   幕墙外环梁\n');
for i = levelPstart1:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1;
    for j = 1:stairColu_num*2	% 每层外筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_f*(i-1); % 
        if j ~= stairColu_num*2
            iN2 = iN1+1;
        else % j = car_num*2 时， 连接的是本环的第一个点，而不是上层内环的第一个点。
            iN2 = iN1+1-stairColu_num*2;
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