%% function
% MGT towerC ramp
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_ramp(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, tube_outerR, levelZaxis, levelPstart, iNO_towerS_init)    %ע�������levelPstart��1x2����
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

% �����µ�ÿ���������б仯������Ҫѭ��(��������)��(��ʱδ���������仯) % ��17600��ʼ���¡����Ͻ���
ramp_inner_R = 12200; % �µ���Ȧ�뾶
ramp_outer_R = 14400; % �µ���Ȧ�뾶

iNO_init = iNO;
levelPstart1 = levelPstart(1); levelPstart2 = levelPstart(2);
lengthlevelZaxis = length(levelZaxis(:));
% XYcoor_i = zeros(car_num,2);	% ��ͲXoY�����1(X)��2(Y)�С�
XYcoor_o = zeros(car_num*2,2);	% ��ͲXoY�����1(X)��2(Y)�С�
XYcoor_ramp_i = zeros(1,car_num*2,2);	% �µ��ڲ�XoY�����1(X)��2(Y)�С�
XYcoor_ramp_o = zeros(1,car_num*2,2);	% �µ����XoY�����1(X)��2(Y)�С�

car_num2pi = 2*pi/car_num;  % speed up

XYcoor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % ����Y��ģ����Ͳһ�� X % Y��ģ�鿴Parts�ļ����ڵ�Ͳ�ļ�
XYcoor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y

XYcoor_o_1(1,1) = sqrt(tube_outerR^2 - XYcoor_i_1(1,2)^2);        % ����Y��ģ����Ͳһ�� X1 ע����Ͳ16���㲢���ǵȽǶȵȷ֡�
XYcoor_o_1(1,2) = XYcoor_i_1(1,2);                                % Y1
XYcoor_o_1(2,:) = coorMir(XYcoor_o_1(1,:), [0,0], XYcoor_i_1);     % X2,Y2
for i = 0:(car_num-1)   % ���������� % ��ת�ֲ��Ƕ�+����Ƕ�
    [XYcoor_o(i*2+1,1), XYcoor_o(i*2+1,2)] = coorTrans(XYcoor_o_1(1,1), XYcoor_o_1(1,2), -car_num2pi*i+Deg_tower); % ��Ͳ������1
    [XYcoor_o(i*2+2,1), XYcoor_o(i*2+2,2)] = coorTrans(XYcoor_o_1(2,1), XYcoor_o_1(2,2), -car_num2pi*i+Deg_tower); % ��Ͳ������2
end

XYcoor_ramp_i_1(1,1) = sqrt(ramp_inner_R^2 - XYcoor_i_1(1,2)^2);        % Y���µ��ڲ���һ�� X1 ע��16���㲢���ǵȽǶȵȷ֡�
XYcoor_ramp_i_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_ramp_i_1(2,:) = coorMir(XYcoor_ramp_i_1(1,:), [0,0], XYcoor_i_1);% ��Y�������߾��񣬵õ�Y��ģ�������֧��X2,Y2
XYcoor_ramp_o_1(1,1) = sqrt(ramp_outer_R^2 - XYcoor_i_1(1,2)^2);        % Y���µ������һ�� X1 ע��16���㲢���ǵȽǶȵȷ֡�
XYcoor_ramp_o_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_ramp_o_1(2,:) = coorMir(XYcoor_ramp_o_1(1,:), [0,0], XYcoor_i_1);% ��Y�������߾��񣬵õ�Y��ģ�������֧��X2,Y2
for j = 1 % �ݶ�1�������µ��������Ⱥ͸߶��йأ����ڻ�����ġ�
    for i = 0:(car_num-1)   % ���������� % ��ת�ֲ��Ƕ�+����Ƕ�
        [XYcoor_ramp_i(j,i*2+1,1), XYcoor_ramp_i(j,i*2+1,2)] = coorTrans(XYcoor_ramp_i_1(1,1), XYcoor_ramp_i_1(1,2), -car_num2pi*i+Deg_tower); % �µ��ڲ������1
        [XYcoor_ramp_i(j,i*2+2,1), XYcoor_ramp_i(j,i*2+2,2)] = coorTrans(XYcoor_ramp_i_1(2,1), XYcoor_ramp_i_1(2,2), -car_num2pi*i+Deg_tower); % �µ��ڲ������2
        [XYcoor_ramp_o(j,i*2+1,1), XYcoor_ramp_o(j,i*2+1,2)] = coorTrans(XYcoor_ramp_o_1(1,1), XYcoor_ramp_o_1(1,2), -car_num2pi*i+Deg_tower); % �µ���������1
        [XYcoor_ramp_o(j,i*2+2,1), XYcoor_ramp_o(j,i*2+2,2)] = coorTrans(XYcoor_ramp_o_1(2,1), XYcoor_ramp_o_1(2,2), -car_num2pi*i+Deg_tower); % �µ���������2
    end
end
% �ֲ�����ϵ ת���� ��������ϵ
% XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_tower(1);
% XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_tower(2);
XYcoor_o(:,1) = XYcoor_o(:,1) + CoC_tower(1);
XYcoor_o(:,2) = XYcoor_o(:,2) + CoC_tower(2);
XYcoor_ramp_i(:,:,1) = XYcoor_ramp_i(:,:,1) + CoC_tower(1);
XYcoor_ramp_i(:,:,2) = XYcoor_ramp_i(:,:,2) + CoC_tower(2);
XYcoor_ramp_o(:,:,1) = XYcoor_ramp_o(:,:,1) + CoC_tower(1);
XYcoor_ramp_o(:,:,2) = XYcoor_ramp_o(:,:,2) + CoC_tower(2);

[~,lengthXYcoor_ramp,~] = size(XYcoor_ramp_i);  % �µ�ÿ��ڵ���
% �µ���������ߵ�ȷ�����ݶ�һȦ�µ�6000��
!theta_temp = arctan ( (XYcoor_ramp_i_1(1,2)+XYcoor_ramp_o_1(1,2)) / (XYcoor_ramp_i_1(1,1)+XYcoor_ramp_o_1(1,1)) );


for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for j = 1:lengthXYcoor_ramp
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_o(j,1),XYcoor_o(j,2),levelZaxis(i));   % ��Ͳ X & Y
    end
    for j = 1:lengthXYcoor_ramp
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_ramp_i(i,j,1),XYcoor_ramp_i(i,j,2),levelZaxis(i));   % ��Ͳ X & Y
    end
    for j = 1:lengthXYcoor_ramp
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_ramp_o(i,j,1),XYcoor_ramp_o(i,j,2),levelZaxis(i));   % ��Ͳ X & Y
    end
end
iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_colu = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% ��Ͳ����iPRO = 1 ������1��
fprintf(fileID,'; Ļǽ��(������Ͳ��)\n');
ELE_iPRO = 1;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart1:(lengthlevelZaxis-1)	% length(A(:)) A����Ԫ�ظ��� % levelPstart �ڼ��㿪ʼͣ�������¼��㿪��
    for j = 1:car_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_ramp*(i-1); % ����Ļǽ���������Ĳ�ͬ����Ļǽ�ڵ�Ϊ������ţ���һ��ֻ��lengthXYcoor_f���ڵ�
        iN2 = iN1+lengthXYcoor_ramp;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) beams
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% ����������iPRO = 3 ������3��
fprintf(fileID,'; ��������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
iNO_towerS = iNO_towerS_init; % ��ʼ��iNO_towerS
for i = levelPstart1:(levelPstart2-1)	% ��ͣ���㣬������Ͳ���
    for j = 1:car_num	% ÿ����Ͳ�Ľڵ���
        for k = 1:2 % һ����Ͳ������������Ͳ������������
            iEL = iEL+1;
            iN1 = iNO_towerS+j+lengthXYcoor2*(i-1); % ����Ϊ��λ������¥�Ľڵ�(��Ͳ)
            iN2 = iNO+lengthXYcoor_ramp*(i-1)+(j-1)*2+k;    % �鵽Ļǽ��Ͳ��0����ٶ�λ�������
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % ����Ԫ�������ڵ��
                ELE_ANGLE, ELE_iSUB);
        end
    end
end
for i = levelPstart2:lengthlevelZaxis	% ͣ���㣬������Ͳ���
    for j = 1:car_num	% ÿ����Ͳ�Ľڵ���
        for k = 1:2 % ������Ͳ������������
            iEL = iEL+1;
            iN1 = iNO_towerS+car_num+lengthXYcoor2*(i-1)+(j-1)*2+k; % ����Ϊ��λ������¥�Ľڵ�(��Ͳ)
            iN2 = iNO+lengthXYcoor_ramp*(i-1)+(j-1)*2+k;    % �鵽Ļǽ��Ͳ��0����ٶ�λ�������
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % ����Ԫ�������ڵ��
                ELE_ANGLE, ELE_iSUB);
        end
    end
end

% ���δ�����iPRO = 4 ������4��
fprintf(fileID,'; ���δ���\n');
ELE_iPRO = 4;
iNO = iNO_init; % ��ʼ��iNO
% �⻷��
fprintf(fileID,';   Ļǽ�⻷��\n');
for i = levelPstart1:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1;
    for j = 1:car_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_ramp*(i-1); % 
        if j ~= car_num*2
            iN2 = iN1+1;
        else % j = car_num*2 ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2 = iN1+1-car_num*2;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
iEL_end = iEL;
fprintf(fileID,'\n');
end