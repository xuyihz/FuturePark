%% function
% MGT tower 2 facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_facade_S2(fileID, iNO, iEL, column_num, CoC_tower, Deg_tower, towerS_column_coor, facade_tower_R, levelZaxis, levelPstart, ~, ~)    %ע�������levelPstart��1x2����
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

levelPstart1 = levelPstart(1);
lengthlevelZaxis = length(levelZaxis(:));
xtemp = towerS_column_coor(1); ytemp = towerS_column_coor(2);
XYcoor_o3 = zeros(lengthlevelZaxis,column_num*2,2);	% ��ͲXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

for j = levelPstart1:lengthlevelZaxis
    XY_o_1x = sqrt(facade_tower_R(j)^2 - ytemp^2);	% Ļǽ��1�� x��
    XY_o_2y = sqrt(facade_tower_R(j)^2 - xtemp^2);  % Ļǽ��2�� y��
    
    XYcoor_o3(j,:,:) = [XY_o_1x, ytemp; xtemp, XY_o_2y; -xtemp, XY_o_2y; -XY_o_1x, ytemp;...
        -XY_o_1x, -ytemp; -xtemp, -XY_o_2y; xtemp, -XY_o_2y; XY_o_1x, -ytemp;];
    if Deg_tower ~= 0
        for i = 1:column_num*2
            XYcoor_o3(j,i,:) = coorTrans(XYcoor_o3(j,i,:), Deg_tower);
        end
    end
end

% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_tower(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_tower(2);

lengthXYcoor_f = column_num*2;  % Ļǽÿ��ڵ���

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for j = 1:lengthXYcoor_f % Ļǽ
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));   % ��Ͳ X & Y
    end
end
iNO_end = iNO;
iEL_end = iEL;
fprintf(fileID,'\n');
end