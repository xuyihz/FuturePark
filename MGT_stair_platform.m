%% function
% MGT stairs
%
% Xu Yi, 25th March 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_stair_platform(fileID, iNO, iEL, CoC_stair, Deg_stair, levelZaxis, stairColu_num, stairL, stairW, stairB, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis = length(levelZaxis(:));

XYcoor_i = zeros(stairColu_num,2);   % ��ͲXoY�����1(X)��2(Y)�С�
XYcoor_pf = XYcoor_i;	% ¥��ƽ̨/��Ϣƽ̨ XoY�����1(X)��2(Y)�С�

% ��Ͳ��
stairXYtemp1 = [stairL/2, stairW/2; -stairL/2, stairW/2]; % ԭʼ���꣬δ��ת��δת����������ϵ
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
% ¥��ƽ̨/��Ϣƽ̨
stairXYtemp1 = [stairL/2+stairB, stairW/2; -stairL/2-stairB, stairW/2]; % ԭʼ���꣬δ��ת��δת����������ϵ
stairXYtemp2 = zeros(length(stairXYtemp1),2);
for i = 1:length(stairXYtemp1)
    stairXYtemp2(i,:) = stairXYtemp1((length(stairXYtemp1)-i+1),:); % ��ʱ����
end
for i = 1:length(stairXYtemp1)
    stairXYtemp2(i,2) = -stairXYtemp1((length(stairXYtemp1)-i+1),2); % ��X��Գ�
end
stairXY_pf = [stairXYtemp1; stairXYtemp2];

for i = 1:stairColu_num   % ����������
    [XYcoor_pf(i,:)] = coorTrans(stairXY_pf(i,:), Deg_stair); % ��Ͳ
end

% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_stair(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_stair(2);
XYcoor_pf(:,1) = XYcoor_pf(:,1) + CoC_stair(1);
XYcoor_pf(:,2) = XYcoor_pf(:,2) + CoC_stair(2);

lengthXYcoor_i = length(XYcoor_i); % ��Ͳÿ��ڵ���
lengthXYcoor_pf = length(XYcoor_pf); % ƽ̨ÿ��ڵ���
lengthXYcoor_all = lengthXYcoor_i + lengthXYcoor_pf;  % ÿ��ڵ�������

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for j = 1:lengthXYcoor_i % �ڲ�4����
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
    end
    for j = 1:lengthXYcoor_pf % �ⲿ4����
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_pf(j,1),XYcoor_pf(j,2),levelZaxis(i));
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
fprintf(fileID,'; ¥������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis	%
    if rem(i,2) ~= 0 % i��������
        iN1 = iNO+1+lengthXYcoor_all*(i-1); % ��Ͳ��1
        iN2 = iN1+3;    % ��Ͳ��4
    else % i��ż����
        iN1 = iNO+2+lengthXYcoor_all*(i-1); % ��Ͳ��2
        iN2 = iN1+1;    % ��Ͳ��3
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % ����Ԫ�������ڵ��
        ELE_ANGLE, ELE_iSUB);
end

fprintf(fileID,'; ¥��������\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis	%
    if rem(i,2) ~= 0 % i��������
        iN_i(1) = iNO+1+lengthXYcoor_all*(i-1); % ��Ͳ��1
        iN_i(2) = iN_i(1)+3;    % ��Ͳ��4
        iN_pf(1) = iNO+lengthXYcoor_i+1+lengthXYcoor_all*(i-1); % ƽ̨��1
        iN_pf(2) = iN_pf(1)+3;    % ƽ̨��4
    else % i��ż����
        iN_i(1) = iNO+2+lengthXYcoor_all*(i-1); % ��Ͳ��2
        iN_i(2) = iN_i(1)+1;    % ��Ͳ��3
        iN_pf(1) = iNO+lengthXYcoor_i+2+lengthXYcoor_all*(i-1); % ƽ̨��2
        iN_pf(2) = iN_pf(1)+1;    % ƽ̨��3
    end
    for j = 1:2 % ����������
        iEL = iEL+1;
        iN1 = iN_i(j);
        iN2 = iN_pf(j);
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end

fprintf(fileID,'; ¥��б��\n');
iNO = iNO_init; % ��ʼ��iNO
iN_i = zeros(1,lengthXYcoor_i);
iN_table = [1,2; 4,3; 2,1; 3,4]; % ������
for i = 1:(lengthlevelZaxis-1)	% б������-1.(�߼�ͬ��)
    for k = 1:lengthXYcoor_i % ��Ͳ4����
        iN_i(k) = iNO+k+lengthXYcoor_all*(i-1); % ��Ͳ��1/2/3/4
    end
    for j = 1:2 % ����б��
    if rem(i,2) ~= 0 % i��������
        iN1 = iN_i( iN_table(j,1) ); % ��Ͳ��1/4
        iN2 = iN_i( iN_table(j,2) ) + lengthXYcoor_all;    % �ϲ���Ͳ��2/3
    else % i��ż����
        iN1 = iN_i( iN_table(j+2,1) ); % ��Ͳ��2/3
        iN2 = iN_i( iN_table(j+2,2) ) + lengthXYcoor_all;    % �ϲ���Ͳ��1/4
    end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end

% ���δ�����iPRO = 4 ������4��
fprintf(fileID,'; ƽ̨�����\n');
ELE_iPRO = 4;
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis	%
    if rem(i,2) ~= 0 % i��������
        iN1 = iNO+lengthXYcoor_i+1+lengthXYcoor_all*(i-1); % ƽ̨��1
        iN2 = iN1+3;    % ƽ̨��4
    else % i��ż����
        iN1 = iNO+lengthXYcoor_i+2+lengthXYcoor_all*(i-1); % ƽ̨��2
        iN2 = iN1+1;    % ƽ̨��3
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
    if rem(i,2) ~= 0    % ������ 1/2+/3+/4
        iN1 = iNO +1 + lengthXYcoor_all*(i-1); % ����1
        iN2 = iN1 +1 + lengthXYcoor_all; % �ϲ�2
        iN3 = iN2 +1; % �ϲ�3
        iN4 = iN1 +3; % ����4
    else % ż���� 1+/2/3/4+
        iN1 = iNO +1 + lengthXYcoor_all*i; % �ϲ�1
        iN2 = iNO +2 + lengthXYcoor_all*(i-1); % ����2
        iN3 = iN2 +1; % ����3
        iN4 = iN1 +3; % �ϲ�4
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2, iN3, iN4,...    % �嵥Ԫ���ĸ��ڵ��
        ELE_iSUB, ELE_iWID);
end

fprintf(fileID,'; 1���¥��ƽ̨��\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis % ƽ̨����Ҫ-1
    if rem(i,2) ~= 0    % ������ 1/4/4'/1'
        iN1 = iNO +1 + lengthXYcoor_all*(i-1); % ��Ͳ1
        iN2 = iN1 +3; % ��Ͳ4
        iN3 = iN2 + lengthXYcoor_i; % ƽ̨4
        iN4 = iN1 + lengthXYcoor_i; % ƽ̨1
    else % ż���� 2/2'/3'/3
        iN1 = iNO +2 + lengthXYcoor_all*(i-1); % ��Ͳ2
        iN2 = iN1 + lengthXYcoor_i; % ƽ̨2
        iN3 = iN2 +1; % ƽ̨3
        iN4 = iN1 +1; % ��Ͳ3
    end
    iEL = iEL+1;
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

fprintf(fileID,'; 1���¥��б�����\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:(lengthlevelZaxis-1) % ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ 1/2+/3+/4
        iN1 = iNO +1 + lengthXYcoor_all*(i-1); % ����1
        iN2 = iN1 +1 + lengthXYcoor_all; % �ϲ�2
        iN3 = iN2 +1; % �ϲ�3
        iN4 = iN1 +3; % ����4
    else % ż���� 1+/2/3/4+
        iN1 = iNO +1 + lengthXYcoor_all*i; % �ϲ�1
        iN2 = iNO +2 + lengthXYcoor_all*(i-1); % ����2
        iN3 = iN2 +1; % ����3
        iN4 = iN1 +3; % �ϲ�4
    end
    fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
end

fprintf(fileID,'; 1���¥��ƽ̨�����\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis % ƽ̨����Ҫ-1
    if rem(i,2) ~= 0    % ������ 1/4/4'/1'
        iN1 = iNO +1 + lengthXYcoor_all*(i-1); % ��Ͳ1
        iN2 = iN1 +3; % ��Ͳ4
        iN3 = iN2 + lengthXYcoor_i; % ƽ̨4
        iN4 = iN1 + lengthXYcoor_i; % ƽ̨1
    else % ż���� 2/2'/3'/3
        iN1 = iNO +2 + lengthXYcoor_all*(i-1); % ��Ͳ2
        iN2 = iN1 + lengthXYcoor_i; % ƽ̨2
        iN3 = iN2 +1; % ƽ̨3
        iN4 = iN1 +1; % ��Ͳ3
    end
    fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

end