%% function
% MGT stairs
%
% Xu Yi, 25th March 2018
% Xu Yi, 28th March 2018, revised

%%
function [iNO_end, iEL_end] = MGT_stair(fileID, iNO, iEL, CoC_tower, levelZaxis, levelPstart, stairN_num, stairL, stairW, CAR, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

stairN_num = 8;
iNO_init = iNO;
XYcor = zeros(stairN_num,4);   % ��ͲXoY�����1(X)��2(Y)��(1~4��)(ע��5~8�о�Ϊ0����ʹ��)�����ƤXoY�����3(X)��4(Y)��(1~8��)��

stairL = 3150; % ¥�ݳ�������̤��ǰ������
stairW = 3250; % ¥�ݿ�
stairXY = [stairL/2, stairW/2; -stairL/2, stairW/2; -stairL/2, -stairW/2; stairL/2, -stairW/2]; % ԭʼ���꣬δ��ת��δת����������ϵ
stairDeg = -pi/6; % ¥����ת�Ƕȣ��ݶ� ˳ʱ��ת30(pi/6)
% ����ϵת�� function coorTrans
% x' = xcos(a) + ysin(a);
% y' = ycos(a) - xsin(a);
for i = 1:stairN_num/2   % ����������
    [XYcor(i,1), XYcor(i,2)] = coorTrans(stairXY(i,1), stairXY(i,2), stairDeg); % ��Ͳ
end
% ��Ͳ����stairXY2 �ݶ�����Ͳ����2000.
m2 = 2000;
stairXY2 = [stairL/2+m2, stairW/2; stairL/2, stairW/2+m2; -stairL/2, stairW/2+m2; -stairL/2-m2, stairW/2;...
    -stairL/2-m2, -stairW/2; -stairL/2, -stairW/2-m2; stairL/2, -stairW/2-m2; stairL/2+m2, -stairW/2];
for i = 1:stairN_num   % ����������
    [XYcor(i,3), XYcor(i,4)] = coorTrans(stairXY2(i,1), stairXY2(i,2), stairDeg); % ��Ͳ
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcor(:,1) = XYcor(:,1) + CoC_tower(1);
XYcor(:,2) = XYcor(:,2) + CoC_tower(2);
XYcor(:,3) = XYcor(:,3) + CoC_tower(1);
XYcor(:,4) = XYcor(:,4) + CoC_tower(2);
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for k = 1:2 % ��Ͳ����Ͳ
        if k == 1 % �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻��ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
            for j = 1:stairN_num/2 % �ڲ�4������
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcor(j,1),XYcor(j,2),levelZaxis(i));
            end
        else
            for j = 1:stairN_num % �ⲿ8������
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcor(j,3),XYcor(j,4),levelZaxis(i));
            end
        end
    end
end
lengthXYcor2 = stairN_num/2+stairN_num;  % ÿ��Ľڵ����������ڲ�4���㣬�ⲿ8���㡣
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
    for j = 1:stairN_num/2	% ÿ����Ͳ�Ľڵ��� stairN_num/2
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
    for j = 1:stairN_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+(stairN_num/2+j)+lengthXYcor2*(i-1); % ��������Ͳ��ͬ������ +stairN_num/2
        iN2 = iN1+lengthXYcor2;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) beams ������ͣ��Ͳ��ͬ���Ҿ�Ϊͬһ����
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% ����������iPRO = 3 ������3��
fprintf(fileID,'; ��������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1
    for j = 1:stairN_num/2	% ÿ����Ͳ�Ľڵ���
        iN1 = iNO+j+lengthXYcor2*(i-1);
        for k = 1:2
            iN2 = iN1+stairN_num/2-(j-1)+(j-1)*2+(k-1); % ��ͣ��Ͳ����������һ��
            iEL = iEL+1;
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
    for j = 1:stairN_num/2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcor2*(i-1);
        if j ~= stairN_num/2
            iN2 = iN1+1;
        else % j = car_num ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������⻷�ĵ�һ���㡣
            iN2 = iN1+1-stairN_num/2;
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
    for j = 1:stairN_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+stairN_num/2+j+lengthXYcor2*(i-1); % �������ڻ�����ͬ�������stairN_num/2
        if j ~= stairN_num
            iN2 = iN1+1;
        else % j = car_num ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2 = iN1+1-stairN_num;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(planner) floor
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_floor = iEL;
ELE_TYPE = 'PLATE'; ELE_iMAT = 2; ELE_iSUB = 2; ELE_iWID = 0; % iMAT = 2���ϻ�����C30 % iSUB = 2 ����

% ���1��iPRO = 2 ������2��
fprintf(fileID,'; 1���¥�ݰ�\n');
ELE_iPRO = 2;
iNO = iNO_init; % ��ʼ��iNO
for i = 2:lengthlevelZaxis % ����ͬ�⻷�� % ������������ͬ��i��ʼΪ2.�����㿪ʼ�� % ��ÿ��ֻ������һ��壬��û��j��ѭ��
    iEL = iEL+1;
    iN1 = iNO+1+lengthXYcor2*(i-1); % ��ʱ��������ĸ���
    iN2 = iN1+1;
    iN3 = iN2+1;
    iN4 = iN3+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2, iN3, iN4,...    % �嵥Ԫ���ĸ��ڵ��
        ELE_iSUB, ELE_iWID);
end
fprintf(fileID,'\n');

%% FLOORLOAD
fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');

LTNAME = ROOF; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0; % ¥�ݺ��� Ŀǰ�������������4/3.5�����0.��ÿ2.1һ��壬�ʿ��ܺ�4.2һ��7/3.5��ࡣ(������)
DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';

iNO = iNO_init; % ��ʼ��iNO
for i = 2:lengthlevelZaxis % ����ͬ�⻷�� % ������������ͬ��i��ʼΪ2.�����㿪ʼ�� % ��ÿ��ֻ������һ��壬��û��j��ѭ��
    iN1 = iNO+1+lengthXYcor2*(i-1); % ��ʱ��������ĸ���
    iN2 = iN1+1;
    iN3 = iN2+1;
    iN4 = iN3+1;
    fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
        LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
        iN1, iN2, iN3, iN4);
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % ��ʼ��iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+stairN_num/2);
CONSTRAINT = 111111; % 6�����ɶ�ȫԼ��
fprintf(fileID,'   %s, %d, \n',...
                NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end