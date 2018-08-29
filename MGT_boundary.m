%% function
% MGT boundary beams.
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_boundary(fileID, iNO, iEL, ~, CoC_tower, ~, ~, ~, levelZaxis, ~, Roof_boundary, ~,...
    CoC_towerS2, CoC_towerS3, CoC_elevator4, CoC_stair5, CoC_stair6, facade_tower2_R, facade_tower3_R, facade_ele4_R, facade_stair5_R, facade_stair6_R,...
    levelTaxis, levelTSaxis, levelSaxis_f)    % 注意这里的levelPstart是1x2数组
%%
CoC = [CoC_tower; CoC_towerS2; CoC_towerS3; CoC_elevator4; CoC_stair5; CoC_stair6];                     % centre of tower1~6
facade_R = {[0,0]; facade_tower2_R; facade_tower3_R; facade_ele4_R; facade_stair5_R; facade_stair6_R};    % facade R of tower1~6
levelZ_f = {levelTaxis; levelTSaxis; levelTSaxis; levelSaxis_f; levelSaxis_f; levelSaxis_f};            % facade Z axis of tower1~6
f_boundary = {[0,0,0,0]; [Roof_boundary(1,:),Roof_boundary(2,:)]; [Roof_boundary(7,:),Roof_boundary(1,:)]; [Roof_boundary(7,:),Roof_boundary(1,:)];...
    [Roof_boundary(2,:),Roof_boundary(3,:),Roof_boundary(3,:),Roof_boundary(4,:)]; [Roof_boundary(4,:),Roof_boundary(5,:),Roof_boundary(5,:),Roof_boundary(6,:)]};  % 1~6塔楼接触的边线
%% NODE 定义商业层和屋面层上坡道投影处的节点
iNO_init1 = iNO;
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
iNO_init2 = iNO;
% 定义每层边线
for i = 2:6 % 目前确定塔楼2~6的边线
    CoC_i = CoC(i,:);                 % centre of tower i
    levelZ_f_i = levelZ_f{i};       % facade Z axis of tower i
    facade_R_i = facade_R{i};       % facade R of tower i
    f_boundary_i = f_boundary{i};   % 塔楼 i 接触的边线
    for j = 1:length(levelZ_f_i) % 层数
        if facade_R_i(j) == 0   % 跳过下面没有幕墙的几层
        else
            boundary_num = length(f_boundary_i)/4; % 一条边线还是两条边线
            for k = 1:boundary_num % 边线数
                k_s = k*4-4;
                f_b_temp1 = [f_boundary_i(k_s+1),f_boundary_i(k_s+2)];
                f_b_temp2 = [f_boundary_i(k_s+3),f_boundary_i(k_s+4)];                
                if coorPerpL(CoC_i, f_b_temp1, f_b_temp2) < facade_R_i(j) % 即边线与圆相交 % coorPerpL是垂线的长度
                    [X_temp, Y_temp, ~] = coorLxCp(CoC_i, facade_R_i(j), f_b_temp1, f_b_temp2); % 相交点
                    for h = 1:2
                        iNO = iNO+1;
                        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：从0度角开始逆时针；从下到上。
                            iNO,X_temp(h),Y_temp(h),levelZ_f_i(j));   % 外筒 X & Y
                    end
                else
                end
            end
        end
    end
end

iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns

%% ELEMENT(frame) beams 定义商业层和屋面层上坡道投影处的环梁
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 屋面环形边梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 屋面环形边梁\n');
ELE_iPRO = 4;
iNO = iNO_init1; % 初始化iNO
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
% 每层环形边梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 每层环形边梁\n');
ELE_iPRO = 4;
iNO = iNO_init2; % 初始化iNO
iN1 = iNO-1;
while iN1 < iNO_end-1
    iEL = iEL+1;
    iN1 = iN1+2;	% 12,34,56...
    iN2 = iN1+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % 梁单元的两个节点号
        ELE_ANGLE, ELE_iSUB);
end
fprintf(fileID,'\n');

iEL_end = iEL;
end
