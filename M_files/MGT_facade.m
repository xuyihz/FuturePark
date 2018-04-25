%% function
% MGT facade
%
% Xu Yi, 28th March 2018
% Xu Yi, 28rd March 2018, revised ȫ��δ��

%%
function [iNO_end, iEL_end] = MGT_facade(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, levelZaxis, levelPstart)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
XYcoor_i = zeros(car_num,2);     % ��ͲXoY�����1(X)��2(Y)�С�
XYcoor_o = zeros(car_num*2,2);	% ��ͲXoY�����1(X)��2(Y)�С�

car_num2pi = 2*pi/car_num;  % speed up

XYcor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % ����Y��ģ����Ͳһ�� X
XYcor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y
XYcor_o_1(1,1) = sqrt(tube_outerR^2 - XYcor_i_1(1,2)^2);        % ����Y��ģ����Ͳһ�� X1 ע����Ͳ16���㲢���ǵȽǶȵȷ֡�
XYcor_o_1(1,2) = XYcor_i_1(1,2);                                % Y1
XYcor_o_1(2,:) = coorMir(XYcor_o_1(1,:), [0,0], XYcor_i_1);     % X2,Y2

for i = 0:(car_num-1)   % ���������� % ��ת�ֲ��Ƕ�+����Ƕ�
    [XYcoor_i(i+1,1), XYcoor_i(i+1,2)] = coorTrans(XYcor_i_1(1), XYcor_i_1(2), -car_num2pi*i+Deg_tower);       % ��Ͳ������
    [XYcoor_o(i*2+1,1), XYcoor_o(i*2+1,2)] = coorTrans(XYcor_o_1(1,1), XYcor_o_1(1,2), -car_num2pi*i+Deg_tower); % ��Ͳ������1
    [XYcoor_o(i*2+2,1), XYcoor_o(i*2+2,2)] = coorTrans(XYcor_o_1(2,1), XYcor_o_1(2,2), -car_num2pi*i+Deg_tower); % ��Ͳ������2
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_tower(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_tower(2);
XYcoor_o(:,1) = XYcoor_o(:,1) + CoC_tower(1);
XYcoor_o(:,2) = XYcoor_o(:,2) + CoC_tower(2);

lengthXYcor2 = length(XYcoor_i(:))/2 + length(XYcoor_o(:))/2;  % ÿ��ڵ���
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for k = 1:2 % ��Ͳ1 ��Ͳ2
        for j = 1:car_num
            iNO = iNO+1;
            if k == 1                                           % �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻��ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
            else
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_o(j*2-1,1),XYcoor_o(j*2-1,2),levelZaxis(i));   % ��Ͳ X1 & Y1
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_o(j*2,1),XYcoor_o(j*2,2),levelZaxis(i));       % ��Ͳ X2 & Y2
            end
        end
    end
end
iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_colu = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% ��Ͳ����iPRO = 2 ������2��
fprintf(fileID,'; ��Ͳ��\n');
ELE_iPRO = 2;
iNO = iNO_init; % ��ʼ��iNO
for i = 1:(lengthlevelZaxis-1)	% length(A(:)) A����Ԫ�ظ���
    for j = 1:car_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcor2*(i-1);
        iN2 = iN1+lengthXYcor2;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end

% ��Ͳ����iPRO = 1 ������1��
fprintf(fileID,'; ��Ͳ��\n');
ELE_iPRO = 1;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart:(lengthlevelZaxis-1)	% length(A(:)) A����Ԫ�ظ��� % levelPstart �ڼ��㿪ʼͣ�������¼��㿪��
    for j = 1:car_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+car_num+j+lengthXYcor2*(i-1); % ��������Ͳ��ͬ������ +car_num
        iN2 = iN1+lengthXYcor2;
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
for i = levelPstart:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1
    for j = 1:car_num	% ÿ����Ͳ�Ľڵ���
        for k = 1:2 % һ����Ͳ������������Ͳ������������
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcor2*(i-1);
            iN2 = iN1-j+car_num+(j-1)*2+k;    % iN1�鵽��Ͳ��0����ټ�car_num�󣬼�Ϊ��ͲY�͵�0��(����Ͳ���һ��)
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
% �ڻ���
fprintf(fileID,';   �ڻ���\n');
for i = 2:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1; ������������ͬ��i��ʼΪ2.�����㿪ʼ�С�
    for j = 1:car_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcor2*(i-1);
        if j ~= car_num
            iN2 = iN1+1;
        else % j = car_num ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������⻷�ĵ�һ���㡣
            iN2 = iN1+1-car_num;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');
% �⻷��
fprintf(fileID,';   �⻷��\n');
for i = levelPstart:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1;
    for j = 1:car_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+car_num+j+lengthXYcor2*(i-1); % �������ڻ�����ͬ�������car_num
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