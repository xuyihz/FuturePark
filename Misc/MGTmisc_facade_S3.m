%% function
% MGT tower 3 facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_facade_S3(fileID, iNO, iEL, column_num, CoC_tower, Deg_tower, towerS_column_coor, facade_tower_R, levelZaxis, levelPstart, ~, ~)    %ע�������levelPstart��1x2����
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

levelPstart1 = levelPstart(1);
lengthlevelZaxis = length(levelZaxis(:));
xtemp = towerS_column_coor(1); ytemp = towerS_column_coor(2);
XYcoor_o3 = zeros(lengthlevelZaxis,column_num*2,2);	% ��ͲXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

Rec_R = sqrt( xtemp^2 + ytemp^2 ); % 4�������γɵľ��ζԽ��ߵ�һ��
for j = levelPstart1:lengthlevelZaxis
    Rtemp = facade_tower_R(j);
    XY_o_2x = xtemp / Rec_R * Rtemp;	% Ļǽ��2�� x��
    XY_o_2y = ytemp / Rec_R * Rtemp;  % Ļǽ��2�� y��
    
    XYcoor_o3(j,:,:) = [Rtemp, 0; XY_o_2x, XY_o_2y; 0, Rtemp; -XY_o_2x, XY_o_2y;...
        -Rtemp, 0; -XY_o_2x, -XY_o_2y; 0, -Rtemp; XY_o_2x, -XY_o_2y;];
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