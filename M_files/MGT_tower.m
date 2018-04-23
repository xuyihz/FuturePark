%% function
% MGT tower
%
% Xu Yi, 19th March 2018
% Xu Yi, 23rd April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_tower(fileID, iNO, iEL, car_num, CoC_tower, tube_innerR, tube_outerR, levelZaxis, levelPstart, CAR, ~, ~)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
XYcor = zeros(car_num,4);   % ��ͲXoY�����1(X)��2(Y)�У���ͲXoY�����3(X)��4(Y)�С�
car_num2pi = 2*pi/car_num;  % speed up
for i = 0:(car_num-1)   % ����������
    XYcor(i+1,1) = tube_innerR * cos(car_num2pi*i);
    XYcor(i+1,2) = tube_innerR * sin(car_num2pi*i);
    XYcor(i+1,3) = tube_outerR * cos(car_num2pi*i);
    XYcor(i+1,4) = tube_outerR * sin(car_num2pi*i);
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcor(:,1) = XYcor(:,1) + CoC_tower(1);
XYcor(:,2) = XYcor(:,2) + CoC_tower(2);
XYcor(:,3) = XYcor(:,3) + CoC_tower(1);
XYcor(:,4) = XYcor(:,4) + CoC_tower(2);
lengthXYcor2 = length(XYcor(:))/2;
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for k = 1:2
        for j = 1:car_num
            iNO = iNO+1;
            if k == 1                                           % �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻��ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcor(j,1),XYcor(j,2),levelZaxis(i));
            else
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcor(j,3),XYcor(j,4),levelZaxis(i));
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
    for j = 1:car_num	% ÿ����Ͳ�Ľڵ���
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
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcor2*(i-1);
        iN2 = iN1+car_num;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
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
    for j = 1:car_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcor2*(i-1)+car_num; % �������ڻ�����ͬ�������car_num
        if j ~= car_num
            iN2 = iN1+1;
        else % j = car_num ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2 = iN1+1-car_num;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) bracings
iEL = bracings(fileID, iNO_init, iEL, car_num, lengthlevelZaxis, levelPstart, lengthXYcor2); % �����Ų��Ǳ�Ҫ���豸���������ݽṹ������Ҫ���������벻�ӡ�

%% ELEMENT(planner) floor
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_floor = iEL;
ELE_TYPE = 'PLATE'; ELE_iMAT = 2; ELE_iSUB = 2; ELE_iWID = 0; % iMAT = 2���ϻ�����C30 % iSUB = 2 ����

% ���1��iPRO = 2 ������2��
fprintf(fileID,'; 1���ͣ����\n');
ELE_iPRO = 2;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart:lengthlevelZaxis % ����ͬ�⻷��
    for j = 1:car_num	% ÿ��ͣ����
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcor2*(i-1); % ��ʱ��������ĸ���
        iN2 = iN1+car_num;
        if j ~= car_num
            iN3 = iN2+1;
            iN4 = iN1+1;
        else % j = car_num ʱ�� ���ӵ������ĵ�һ���㣬�������ϲ�ĵ�һ���㡣
            iN3 = iN2+1-car_num;
            iN4 = iN1+1-car_num;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2, iN3, iN4,...    % �嵥Ԫ���ĸ��ڵ��
            ELE_iSUB, ELE_iWID);
    end
end
iEL_end = iEL;
fprintf(fileID,'\n');

%% FLOORLOAD
fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');

LTNAME = CAR; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0;
DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';

iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart:lengthlevelZaxis % ����ͬ1���ͣ���壬��ͬ�⻷��
    for j = 1:car_num	% ÿ��ͣ����
        iN1 = iNO+j+lengthXYcor2*(i-1); % ��ʱ��������ĸ���
        iN2 = iN1+car_num;
        if j ~= car_num
            iN3 = iN2+1;
            iN4 = iN1+1;
        else % j = car_num ʱ�� ���ӵ������ĵ�һ���㣬�������ϲ�ĵ�һ���㡣
            iN3 = iN2+1-car_num;
            iN4 = iN1+1-car_num;
        end
        fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
    end
end
fprintf(fileID,'\n');

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % ��ʼ��iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+car_num);
CONSTRAINT = 111111; % 6�����ɶ�ȫԼ��
fprintf(fileID,'   %s, %d, \n',...
                NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end

%% ELEMENT(frame) bracings
function  iEL = bracings(fileID, iNO_init, iEL, car_num, lengthlevelZaxis, levelPstart, lengthXYcor2)
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_bracing = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% б�ţ�iPRO = 5 ������5��
fprintf(fileID,'; ������\n');
ELE_iPRO = 5;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart:2:(lengthlevelZaxis-2)	% ����������Ϊÿ����һ������Ϊ���2�� �˴�����������ͬ��������һ�ţ���Ҫ-2
    for j = 1:car_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+(j+car_num)+lengthXYcor2*(i-1); % ����������Ԫ��ͬ
        if j ~= car_num
            iN2 = iN1+1+lengthXYcor2*2;
        else % j = car_num ʱ�� ���ӵ�������⻷�ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2 = iN1+1+lengthXYcor2*2-car_num;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');
end