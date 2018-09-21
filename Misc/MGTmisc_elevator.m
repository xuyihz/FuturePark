%% function
% MGT elevator
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGTmisc_elevator(fileID, iNO, iEL, CoC_elevator, Deg_elevator, facade_ele_R, levelZaxis_f, levelPstart1, elevatorColu_num, elevatorR, ~, ~, ~, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

lengthlevelZaxis_f = length(levelZaxis_f(:));

XYcoor_i = zeros(elevatorColu_num,2);   % 8����(�������м䲻������1�����Բ�ĵ�)��ͲXoY�����1(X)��2(Y)�С�
elevatorColu_o_num = elevatorColu_num; % ��Ͳ8����,��Ͳ8����
XYcoor_o3 = zeros(lengthlevelZaxis_f,elevatorColu_o_num,2);	% 12���㣬��ͲXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

% ��Ͳ��
elevatorXYtemp = elevatorR/sqrt(2); % ������
elevatorXY = [elevatorXYtemp,elevatorXYtemp; -elevatorXYtemp,elevatorXYtemp;...
    -elevatorXYtemp,-elevatorXYtemp; elevatorXYtemp,-elevatorXYtemp;...
    elevatorXYtemp,350; -elevatorXYtemp,350; 0,elevatorXYtemp; 0,350]; % ԭʼ���꣬δ��ת��δת����������ϵ
for i = 1:elevatorColu_num   % ����������
    XYcoor_i(i,:) = coorTrans(elevatorXY(i,:), Deg_elevator); % ��Ͳ
end
% ��ͲelevatorXY2 % ��Ͳ��Ҫ����Ļǽ���Ƥ���߶�λ % �����ظ߶ȵ�ѭ��(��������)
for j = levelPstart1:lengthlevelZaxis_f
    elevatorXY_o_temp = sqrt(facade_ele_R(j)^2 - elevatorXYtemp^2); %
    elevatorXY_o = [elevatorXY_o_temp,elevatorXYtemp; elevatorXYtemp,elevatorXY_o_temp;...
        -elevatorXYtemp,elevatorXY_o_temp; -elevatorXY_o_temp,elevatorXYtemp;...
        -elevatorXY_o_temp,-elevatorXYtemp; -elevatorXYtemp,-elevatorXY_o_temp;...
        elevatorXYtemp,-elevatorXY_o_temp; elevatorXY_o_temp,-elevatorXYtemp];
    
    for i = 1:elevatorColu_o_num   % ����������
        [XYcoor_o3(j,i,:)] = coorTrans(elevatorXY_o(i,:), Deg_elevator); % ��Ͳ
    end
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_elevator(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_elevator(2);

for i = 1:lengthlevelZaxis_f  % length(A(:)) A����Ԫ�ظ���
    % ��Ͳ����Ͳ % �ڵ��Ź�����ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
    for j = 1:elevatorColu_o_num % �ⲿ8������
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis_f(i));
    end
end
iNO_end = iNO;
iEL_end = iEL;
fprintf(fileID,'\n');
end