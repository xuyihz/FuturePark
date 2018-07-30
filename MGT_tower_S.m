%% function
% MGT tower
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_tower_S(fileID, iNO, iEL, column_num, CoC_tower, Deg_tower, towerS_column_coor, levelZaxis, levelPstart1o2, CAR, ~, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
xtemp = towerS_column_coor(1); ytemp = towerS_column_coor(2);
XYcoor_i = [ xtemp,ytemp; -xtemp,ytemp; -xtemp,-ytemp; xtemp,-ytemp ]; % 小塔柱坐标

if Deg_tower == 0 % 旋转角度
else
    for i = 1:column_num
        XYcoor_i(i,:) = coorTrans(XYcoor_i(i,:), Deg_tower);
    end
end

% 局部坐标系 转换至 整体坐标系
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_tower(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_tower(2);

lengthXYcoor2 = length(XYcoor_i(:))/2;  % 每层节点数
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for j = 1:column_num % 节点编号规则：从0度角开始逆时针；从下到上。
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
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
fprintf(fileID,'; 小塔柱\n');
ELE_iPRO = 2;
iNO = iNO_init; % 初始化iNO
for i = 1:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数
    for j = 1:column_num	% 每层内筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor2*(i-1);
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

% 环形次梁；iPRO = 4 截面编号4。
fprintf(fileID,'; 环形次梁\n');
ELE_iPRO = 4;
iNO = iNO_init; % 初始化iNO
% 内环梁
fprintf(fileID,';   内环梁\n');
for i = 2:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1; 此行与主梁不同，i起始为2.即二层开始有。
    for j = 1:column_num	% 每层内筒的节点数
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor2*(i-1);
        if j ~= column_num
            iN2 = iN1+1;
        else % j = car_num 时， 连接的是本环的第一个点，而不是外环的第一个点。
            iN2 = iN1+1-column_num;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 梁单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) bracings 人字撑
% iEL = bracings(fileID, iNO_init, iEL, car_num, lengthlevelZaxis, levelPstart, lengthXYcoor2); % 螺旋撑不是必要的设备构件，根据结构计算需要，决定加与不加。

% %% ELEMENT(planner) floor
% fprintf(fileID,'*ELEMENT    ; Elements\n');
% fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');
% 
% % iEL_init_floor = iEL;
% ELE_TYPE = 'PLATE'; ELE_iMAT = 2; ELE_iSUB = 2; ELE_iWID = 0; % iMAT = 2材料混凝土C30 % iSUB = 2 薄板
% 
% % 板厚1；iPRO = 2 截面编号2。
% fprintf(fileID,'; 1厚板停车板\n');
% ELE_iPRO = 2;
% iNO = iNO_init; % 初始化iNO
% for i = levelPstart1o2:(lengthlevelZaxis-1) % 此行同外环梁，屋面层先不建即(lengthlevelZaxis-1)
%     for j = 1:column_num	% 每层停车数
%         iEL = iEL+1;
%         iN1 = iNO+j+lengthXYcoor2*(i-1);     % 逆时针板四周四个点
%         iN2 = iN1-j+1+column_num+(j-1)*2+1;	% iN1归到内筒第一点后再加car_num后，即为外筒第一点(即Y型第一点，实际起点应再+1，即为Y型第二点)
%         if j ~= column_num
%             iN3 = iN2+1;
%             iN4 = iN1+1;
%         else % j = car_num 时， 连接的是这层的第一个点，而不是上层的第一个点。
%             iN3 = iN2+1-column_num*2;
%             iN4 = iN1+1-column_num;
%         end
%         fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
%             iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
%             iN1, iN2, iN3, iN4,...    % 板单元的四个节点号
%             ELE_iSUB, ELE_iWID);
%     end
% end
iEL_end = iEL;
% fprintf(fileID,'\n');
% 
% %% FLOORLOAD
% fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
% fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');
% 
% LTNAME = CAR; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0;
% DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';
% 
% iNO = iNO_init; % 初始化iNO
% for i = levelPstart1o2:(lengthlevelZaxis-1) % 此行同1厚板停车板，即同外环梁，屋面层先不建即(lengthlevelZaxis-1)
%     for j = 1:column_num	% 每层停车数
%         iN1 = iNO+j+lengthXYcoor2*(i-1); % 逆时针板四周四个点
%         iN2 = iN1-j+1+column_num+(j-1)*2+1;	% iN1归到内筒第一点后再加car_num后，即为外筒第一点(即Y型第一点，实际起点应再+1，即为Y型第二点)
%         if j ~= column_num
%             iN3 = iN2+1;
%             iN4 = iN1+1;
%         else % j = car_num 时， 连接的是这层的第一个点，而不是上层的第一个点。
%             iN3 = iN2+1-column_num*2;
%             iN4 = iN1+1-column_num;
%         end
%         fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
%             LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
%             iN1, iN2, iN3, iN4);
%     end
% end
% fprintf(fileID,'\n');

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % 初始化iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+column_num);
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