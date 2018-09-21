%% function
% MGT stairs
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_stair(fileID, iNO, iEL, CoC_stair, Deg_stair, facade_stair_R, levelZaxis_f, levelPstart1, stairColu_num, stairL, stairW, ~, ~, ~, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

lengthlevelZaxis_f = length(levelZaxis_f(:));

XYcoor_i = zeros(stairColu_num,2);   % ��ͲXoY�����1(X)��2(Y)�С�
XYcoor_o3 = zeros(lengthlevelZaxis_f,stairColu_num*3,2);	% ����*3����ͲXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

% ��Ͳ��
stairXYtemp1 = [stairL/2, stairW/2; -stairL/2, stairW/2]; % ԭʼ���꣬δ��ת��δת����������ϵ
stairXYtemp1Deg = coorDeg([0,0], [1,0], stairXYtemp1(1,:)); % ��Ͳ��һ����X��ĽǶ�
stairXYtemp2 = zeros(length(stairXYtemp1),2);
for i = 1:length(stairXYtemp1)
    stairXYtemp2(i,:) = stairXYtemp1((length(stairXYtemp1)-i+1),:); % ��ʱ����
end
for i = 1:length(stairXYtemp1)
    stairXYtemp2(i,2) = -stairXYtemp1((length(stairXYtemp1)-i+1),2); % ��X��Գ�
end
stairXY_i = [stairXYtemp1; stairXYtemp2];

for i = 1:stairColu_num   % ����������
    [XYcoor_i(i,:)] = coorTrans(stairXY_i(i,:), Deg_stair); % ��Ͳ
end
% Ļǽ��Ͳ
for j = levelPstart1:lengthlevelZaxis_f
    XYcoor_o_x = sqrt(facade_stair_R(j)^2 - (stairW/2)^2);  % Ļǽ��Ȧ�㵱y=��stairW/2ʱ��x�ľ���ֵ
    XYcoor_o_y = sqrt(facade_stair_R(j)^2 - (stairL/2)^2);  % Ļǽ��Ȧ�㵱x=��stairL/2ʱ��y�ľ���ֵ
    
    stairXYtemp = coorTrans([facade_stair_R(j), 0], -stairXYtemp1Deg); % б���˵����ꡣ
    
    stairXY_o_temp1 = [XYcoor_o_x, stairW/2; stairXYtemp; stairL/2, XYcoor_o_y];  % 3�����ģ��
    stairXY_o_temp2 = zeros(length(stairXY_o_temp1),2);
    for i = 1:length(stairXY_o_temp1)
        stairXY_o_temp2(i,:) = stairXY_o_temp1((length(stairXY_o_temp1)-i+1),:); % ��ʱ����
    end
    for i = 1:length(stairXY_o_temp1)
        stairXY_o_temp2(i,1) = -stairXY_o_temp1((length(stairXY_o_temp1)-i+1),1); % ��Y��Գ�
    end
    stairXY_o_temp12 = [stairXY_o_temp1; stairXY_o_temp2];
    stairXY_o_temp4 = zeros(length(stairXY_o_temp12),2);
    for i = 1:length(stairXY_o_temp12)
        stairXY_o_temp4(i,:) = stairXY_o_temp12((length(stairXY_o_temp12)-i+1),:); % ��ʱ����
    end
    for i = 1:length(stairXY_o_temp12)
        stairXY_o_temp4(i,2) = -stairXY_o_temp12((length(stairXY_o_temp12)-i+1),2); % ��X��Գ�
    end
    stairXY_o = [stairXY_o_temp1; stairXY_o_temp2; stairXY_o_temp4]; % Ļǽ��Ͳ��������
    
    for i = 1:length(stairXY_o)   % ����������
        [XYcoor_o3(j,i,:)] = coorTrans(stairXY_o(i,:), Deg_stair); % Ļǽ��Ͳ������ % ����ת����Ƕ�
    end
end

% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_stair(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_stair(2);

for i = 1:lengthlevelZaxis_f  % length(A(:)) A����Ԫ�ظ���
    % ��Ͳ����Ͳ % �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻��ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
    for j = 1:stairColu_num*3 % �ⲿ12������
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis_f(i));
    end
end
iNO_end = iNO;
iEL_end = iEL;
fprintf(fileID,'\n');
end