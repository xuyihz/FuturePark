%% function
% MGT elevator
%
% Xu Yi, 24th April 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_elevator_platform(fileID, iNO, iEL, CoC_elevator, Deg_elevator, levelZaxis, elevatorColu_num, elevatorR, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis = length(levelZaxis(:));

XYcoor_i = zeros(elevatorColu_num,2);   % 8����(�������м䲻������1�����Բ�ĵ�)��ͲXoY�����1(X)��2(Y)�С�

% ��Ͳ��
elevatorXYtemp = elevatorR/sqrt(2);
elevatorXY = [elevatorXYtemp,elevatorXYtemp; -elevatorXYtemp,elevatorXYtemp;...
    -elevatorXYtemp,-elevatorXYtemp; elevatorXYtemp,-elevatorXYtemp;...
    elevatorXYtemp,350; -elevatorXYtemp,350; 0,elevatorXYtemp; 0,350]; % ԭʼ���꣬δ��ת��δת����������ϵ
for i = 1:elevatorColu_num   % ����������
    XYcoor_i(i,:) = coorTrans(elevatorXY(i,:), Deg_elevator); % ��Ͳ
end

% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_elevator(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_elevator(2);

lengthXYcoor_i = length(XYcoor_i); % ��Ͳÿ��ڵ���
lengthXYcoor_all = lengthXYcoor_i;  % ÿ��ڵ�������

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    % ��Ͳ ���µ��ϡ�
    for j = 1:lengthXYcoor_i % �ڲ�7������(��һ��������)
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
    end
end
iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) beams ������ͣ��Ͳ��ͬ���Ҿ�Ϊͬһ����
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% ������iPRO = 3 ������3��
fprintf(fileID,'; ����Ͳ¥��б��\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
iN_i = zeros(1,lengthXYcoor_i);
iN_table = [5,6; 4,3; 6,5; 3,4]; % ������ % ÿ��Ϊ�����˽ڵ��
for i = 1:(lengthlevelZaxis-1)	% б������-1.(�߼�ͬ��)
    for k = 1:lengthXYcoor_i % ��Ͳ8����
        iN_i(k) = iNO+k+lengthXYcoor_all*(i-1); % ��Ͳ��1~8
    end
    for j = 1:2 % ����б��
        if rem(i,2) ~= 0 % i��������
            iN1 = iN_i( iN_table(j,1) ); % ��Ͳ��5/4
            iN2 = iN_i( iN_table(j,2) ) + lengthXYcoor_all;    % �ϲ���Ͳ��6/3
        else % i��ż����
            iN1 = iN_i( iN_table(j+2,1) ); % ��Ͳ��6/3
            iN2 = iN_i( iN_table(j+2,2) ) + lengthXYcoor_all;    % �ϲ���Ͳ��5/4
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end

fprintf(fileID,'; ����Ͳ¥�ݷ����\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis	%
    if rem(i,2) ~= 0 % i��������
        iN1 = iNO+4+lengthXYcoor_all*(i-1); % ��Ͳ��4
        iN2 = iN1+1;    % ��Ͳ��5
    else % i��ż����
        iN1 = iNO+3+lengthXYcoor_all*(i-1); % ��Ͳ��3
        iN2 = iN1+3;    % ��Ͳ��6
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % ����Ԫ�������ڵ��
        ELE_ANGLE, ELE_iSUB);
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
fprintf(fileID,'; 1���¥��б��\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:(lengthlevelZaxis-1) % ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ 5/6+/3+/4
        iN1 = iNO +5 + lengthXYcoor_all*(i-1); % ����5
        iN2 = iN1 +1 + lengthXYcoor_all; % �ϲ�6
        iN3 = iN2 -3; % �ϲ�3
        iN4 = iN1 -1; % ����4
    else % ż���� 5+/6/3/4+
        iN1 = iNO +5 + lengthXYcoor_all*i; % �ϲ�5
        iN2 = iNO +6 + lengthXYcoor_all*(i-1); % ����6
        iN3 = iN2 -3; % ����3
        iN4 = iN1 -1; % �ϲ�4
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2, iN3, iN4,...    % �嵥Ԫ���ĸ��ڵ��
        ELE_iSUB, ELE_iWID);
end
fprintf(fileID,'\n');

%% FLOORLOAD ¥����ؼӲ���
fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');

LTNAME = ROOF; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0; % ¥�ݺ��� Ŀǰ�������������4/3.5�����0.��ÿ2.1һ��壬�ʿ��ܺ�4.2һ��7/3.5��ࡣ(������)
DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';

fprintf(fileID,'; 1���¥��б�����\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:(lengthlevelZaxis-1) % ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ 5/6+/3+/4
        iN1 = iNO +5 + lengthXYcoor_all*(i-1); % ����5
        iN2 = iN1 +1 + lengthXYcoor_all; % �ϲ�6
        iN3 = iN2 -3; % �ϲ�3
        iN4 = iN1 -1; % ����4
    else % ż���� 5+/6/3/4+
        iN1 = iNO +5 + lengthXYcoor_all*i; % �ϲ�5
        iN2 = iNO +6 + lengthXYcoor_all*(i-1); % ����6
        iN3 = iN2 -3; % ����3
        iN4 = iN1 -1; % �ϲ�4
    end
    fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

end