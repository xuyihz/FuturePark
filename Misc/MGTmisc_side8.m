%% function
% MGT side tower 8
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_side8(fileID, iNO, iEL, CoC_side, facade_side_R, levelZaxis, levelPstart, Roof_boundary, ~, ~, ~, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

sideRadius = 5250;
Corner_coor = [Roof_boundary(6,:);Roof_boundary(7,:)]; % �����ǵ�

[X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor(1,:), Roof_boundary(5,:)); %����ߵ�һ�� % [X, Y, Len] = coorLxCp(C0, R, P1, P2);
if X_temp(1) < X_temp(2)
    XYcoor_1 = [X_temp(1),Y_temp(1)];
else
    XYcoor_1 = [X_temp(2),Y_temp(2)];
end
[X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor(2,:), Roof_boundary(1,:)); %�ұߵ�һ�� % [X, Y, Len] = coorLxCp(C0, R, P1, P2);
if X_temp(1) < X_temp(2)
    XYcoor_3 = [X_temp(1),Y_temp(1)];
else
    XYcoor_3 = [X_temp(2),Y_temp(2)];
end
[X_temp, Y_temp, ~] = coorLxCp(CoC_side, sideRadius, Corner_coor(1,:), CoC_side);
if X_temp(1) < X_temp(2)
    XYcoor_2 = [X_temp(1),Y_temp(1)];
else
    XYcoor_2 = [X_temp(2),Y_temp(2)];
end

XYcoor_i = [XYcoor_1; XYcoor_2; XYcoor_3]; % ��ͲXoY�����1(X)��2(Y)�С� % ���ǵ�
sideColu_num = length(XYcoor_i);

lengthlevelZaxis = length(levelZaxis(:));
XYcoor_o3 = zeros(lengthlevelZaxis,sideColu_num,2); % ���ƤXoY�����1(X)��2(Y)�С�
for i = 1:(levelPstart-1)
    for j = 1:sideColu_num
        XYcoor_o3(i,j,:) = CoC_side;
    end
end
for i = levelPstart:lengthlevelZaxis % levelPstart���¶�Ϊ0
    for j = 1:sideColu_num
        if j == sideColu_num
            k = 2;
        else
            k = 1;
        end
        [X_temp, Y_temp, ~] = coorLxCp(CoC_side, facade_side_R(i), Corner_coor(k,:), XYcoor_i(j,:));
        if X_temp(1) < X_temp(2)
            XYcoor_temp = [X_temp(1),Y_temp(1)];
        else
            XYcoor_temp = [X_temp(2),Y_temp(2)];
        end
        XYcoor_o3(i,j,:) = XYcoor_temp;
    end
end

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for j = 1:sideColu_num % �ⲿ����
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));
    end
end
iNO_end = iNO;
iEL_end = iEL;
fprintf(fileID,'\n');
end