%% function
% MGT side tower
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_side(fileID, iNO, iEL, CoC_side, facade_side_R, levelZaxis, levelPstart, Roof_boundary, ~, ~, ~, tower_num, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

switch tower_num % 塔标号
    case 7
        sideRadius = 3300;
        Corner_coor = Roof_boundary(4,:);
        
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Roof_boundary(3,:)); %左边第一点 % [X, Y, Len] = coorLxCp(C0, R, P1, P2);
        if X_temp(1) < X_temp(2)
            XYcoor_1 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_1 = [X_temp(2),Y_temp(2)];
        end
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Roof_boundary(5,:)); %右边第一点 % [X, Y, Len] = coorLxCp(C0, R, P1, P2);
        if X_temp(1) > X_temp(2)
            XYcoor_4 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_4 = [X_temp(2),Y_temp(2)];
        end
        Deg = coorDeg(Corner_coor, XYcoor_1, XYcoor_4);
        Temp_coor_trans = coorTransLoc(Corner_coor, XYcoor_1, -Deg/3);
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Temp_coor_trans);
        if X_temp(1) < X_temp(2)
            XYcoor_2 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_2 = [X_temp(2),Y_temp(2)];
        end
        Temp_coor_trans = coorTransLoc(Corner_coor, XYcoor_1, -Deg/3*2);
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Temp_coor_trans);
        if X_temp(1) < X_temp(2)
            XYcoor_3 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_3 = [X_temp(2),Y_temp(2)];
        end
        
        XYcoor_i = [XYcoor_1; XYcoor_2; XYcoor_3; XYcoor_4]; % 内筒XoY坐标第1(X)、2(Y)列。 % 除角点
        sideColu_num = length(XYcoor_i);
        
        lengthlevelZaxis = length(levelZaxis(:));
        XYcoor_o3 = zeros(lengthlevelZaxis,sideColu_num,2); % 外表皮XoY坐标第1(X)、2(Y)列。
        for i = 1:(levelPstart-1)
            for j = 1:sideColu_num
                XYcoor_o3(i,j,:) = CoC_side;
            end
        end
        for i = levelPstart:lengthlevelZaxis % levelPstart以下都为0
            for j = 1:sideColu_num
                [X_temp, Y_temp, ~] = coorLxCp(CoC_side, facade_side_R(i), Corner_coor, XYcoor_i(j,:));
                if Y_temp(1) < Y_temp(2)
                    XYcoor_temp = [X_temp(1),Y_temp(1)];
                else
                    XYcoor_temp = [X_temp(2),Y_temp(2)];
                end
                XYcoor_o3(i,j,:) = XYcoor_temp;
            end
        end
    case 9
        sideRadius = 3300;
        Corner_coor = Roof_boundary(1,:);
        
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Roof_boundary(7,:)); %左边第一点 % [X, Y, Len] = coorLxCp(C0, R, P1, P2);
        if X_temp(1) > X_temp(2)
            XYcoor_1 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_1 = [X_temp(2),Y_temp(2)];
        end
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Roof_boundary(2,:)); %右边第一点 % [X, Y, Len] = coorLxCp(C0, R, P1, P2);
        if X_temp(1) < X_temp(2)
            XYcoor_4 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_4 = [X_temp(2),Y_temp(2)];
        end
        %         Deg = coorDeg(Corner_coor, XYcoor_1, XYcoor_4);
        Temp_coor_trans = coorTransLoc(Corner_coor, XYcoor_1, -pi/4); %该角塔分隔45°，35°，45°
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Temp_coor_trans);
        if X_temp(1) > X_temp(2)
            XYcoor_2 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_2 = [X_temp(2),Y_temp(2)];
        end
        Temp_coor_trans = coorTransLoc(Corner_coor, XYcoor_4, pi/4); %该角塔分隔45°，35°，45°
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Temp_coor_trans);
        if X_temp(1) > X_temp(2)
            XYcoor_3 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_3 = [X_temp(2),Y_temp(2)];
        end
        
        XYcoor_i = [XYcoor_1; XYcoor_2; XYcoor_3; XYcoor_4]; % 内筒XoY坐标第1(X)、2(Y)列。 % 除角点
        sideColu_num = length(XYcoor_i);
        
        lengthlevelZaxis = length(levelZaxis(:));
        XYcoor_o3 = zeros(lengthlevelZaxis,sideColu_num,2); % 外表皮XoY坐标第1(X)、2(Y)列。
        for i = 1:(levelPstart-1)
            for j = 1:sideColu_num
                XYcoor_o3(i,j,:) = CoC_side;
            end
        end
        for i = levelPstart:lengthlevelZaxis % levelPstart以下都为0
            for j = 1:sideColu_num
                [X_temp, Y_temp, ~] = coorLxCp(CoC_side, facade_side_R(i), Corner_coor, XYcoor_i(j,:));
                if j == 1
                    if X_temp(1) > X_temp(2)
                        XYcoor_temp = [X_temp(1),Y_temp(1)];
                    else
                        XYcoor_temp = [X_temp(2),Y_temp(2)];
                    end
                else
                    if Y_temp(1) > Y_temp(2)
                        XYcoor_temp = [X_temp(1),Y_temp(1)];
                    else
                        XYcoor_temp = [X_temp(2),Y_temp(2)];
                    end
                end
                XYcoor_o3(i,j,:) = XYcoor_temp;
            end
        end
    case 10
        sideRadius = 2500;
        Corner_coor = Roof_boundary(2,:);
        
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Roof_boundary(1,:)); %左边第一点 % [X, Y, Len] = coorLxCp(C0, R, P1, P2);
        if X_temp(1) > X_temp(2)
            XYcoor_1 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_1 = [X_temp(2),Y_temp(2)];
        end
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Roof_boundary(3,:)); %右边第一点 % [X, Y, Len] = coorLxCp(C0, R, P1, P2);
        if X_temp(1) > X_temp(2)
            XYcoor_3 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_3 = [X_temp(2),Y_temp(2)];
        end
        Deg = coorDeg(Corner_coor, XYcoor_1, XYcoor_3);
        Temp_coor_trans = coorTransLoc(Corner_coor, XYcoor_1, -Deg/2);
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor, Temp_coor_trans);
        if X_temp(1) > X_temp(2)
            XYcoor_2 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_2 = [X_temp(2),Y_temp(2)];
        end
        
        XYcoor_i = [XYcoor_1; XYcoor_2; XYcoor_3]; % 内筒XoY坐标第1(X)、2(Y)列。 % 除角点
        sideColu_num = length(XYcoor_i);
        
        lengthlevelZaxis = length(levelZaxis(:));
        XYcoor_o3 = zeros(lengthlevelZaxis,sideColu_num,2); % 外表皮XoY坐标第1(X)、2(Y)列。
        for i = 1:(levelPstart-1)
            for j = 1:sideColu_num
                XYcoor_o3(i,j,:) = CoC_side;
            end
        end
        for i = levelPstart:lengthlevelZaxis % levelPstart以下都为0
            for j = 1:sideColu_num
                [X_temp, Y_temp, ~] = coorLxCp(CoC_side, facade_side_R(i), Corner_coor, XYcoor_i(j,:));
                if X_temp(1) > X_temp(2)
                    XYcoor_temp = [X_temp(1),Y_temp(1)];
                else
                    XYcoor_temp = [X_temp(2),Y_temp(2)];
                end
                XYcoor_o3(i,j,:) = XYcoor_temp;
            end
        end
end

for i = 1:lengthlevelZaxis  % length(A(:)) A向量元素个数
    for j = 1:sideColu_num % 外部柱子
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));
    end
end
iNO_end = iNO;
iEL_end = iEL;
fprintf(fileID,'\n');
end