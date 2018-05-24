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
    for k = 1:2 % ��Ͳ��ƽ̨
        if k == 1 % �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻��ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
            for j = 1:lengthXYcoor_i % �ڲ�4����
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
            end
        else
            for j = 1:lengthXYcoor_pf % �ⲿ4����
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_pf(j,1),XYcoor_pf(j,2),levelZaxis(i));
            end
        end
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
        iN_i1 = iNO+1+lengthXYcoor_all*(i-1); % ��Ͳ��1
        iN_i4 = iN_i1+lengthXYcoor_i;    % ��Ͳ��4
        iN_pf1 = ;
        iN_pf4 = ;
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

% ���δ�����iPRO = 4 ������4��
fprintf(fileID,'; ���δ���\n');
ELE_iPRO = 4;
iNO = iNO_init; % ��ʼ��iNO
% �⻷��
fprintf(fileID,';   Ļǽ�⻷��\n');
for i = levelPstart1:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1;
    for j = 1:lengthXYcoor_pf	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+lengthXYcoor_i+j+lengthXYcoor_all*(i-1); % 
        if j ~= slengthXYcoor_f
            iN2 = iN1+1;
        else % j = lengthXYcoor_f ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2 = iN1+1-lengthXYcoor_pf;
        end
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) beams ¥��¥��ƽ̨����Ϣƽ̨
% ������iPRO = 3 ������3��
fprintf(fileID,'; ¥������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO








for i = 1:(lengthlevelZaxis-1)	% ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬������б�����
        iNcon1 = iNO+1+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor_all+1;
        iNcon4 = iNcon1+lengthXYcoor_all+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor_all-1;
        iNcon4 = iNcon1+lengthXYcoor_all+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+stairColu_num*2-1;
    end
    if i < levelPstart1 % ���ǵײ���Ļǽ
        k_end = 2;
    else
        k_end = 4;
    end
    for k = 1:k_end % �������������� % б��+ƽ��
        if k == 1
            iN1 = iNcon1; iN2 = iNcon3;
        elseif k == 2
            iN1 = iNcon2; iN2 = iNcon4;
        elseif k == 3
            iN1 = iNcon3; iN2 = iNcon5;
        elseif k == 4
            iN1 = iNcon4; iN2 = iNcon6;
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'; ¥�ݿ�������\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1 % ÿ��һ���ᴩ��������
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬��������Ͳ�������
        iNcon1 = iNO+1+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+3;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+1;
    end
    if i < levelPstart1 % ���ǵײ���Ļǽ
        k_end = 1;
    else
        k_end = 3;
    end
    for k = 1:k_end % ��������
        if k == 1
            iN1 = iNcon1;
            iN2 = iNcon2;
        elseif k == 2
            iN1 = iNcon1;
            iN2 = iN1 + 5;
        elseif k == 3
            iN1 = iNcon2;
            iN2 = iN1 + 7;
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end

% ��������iPRO = 4 ������4��
fprintf(fileID,'; ���δ���\n');
ELE_iPRO = 4;
iNO = iNO_init; % ��ʼ��iNO
% �⻷�� % �ο�¥�ݳ�������
fprintf(fileID,';   �⻷��\n');
for i = levelPstart1:(lengthlevelZaxis-1)	% ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬������б�����
        iNcon1 = iNO+6+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+5;
        iNcon3 = iNcon1+lengthXYcoor_all+1;
        iNcon4 = iNcon2+lengthXYcoor_all-1;
        iNcon5 = iNcon3+1;
        iNcon6 = iNcon4-1;
    else % ż����
        iNcon1 = iNO+7+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor_all-1;
        iNcon4 = iNcon2+lengthXYcoor_all+1;
        iNcon5 = iNcon3-1;
        iNcon6 = iNcon4+1;
    end
    for k = 1:5 % ���г������������Σ�����һ��һ�� % б��+ƽ��
        if k == 1
            iN1 = iNcon1; iN2 = iNcon3;
        elseif k == 2
            iN1 = iNcon2; iN2 = iNcon4;
        elseif k == 3
            iN1 = iNcon3; iN2 = iNcon5;
        elseif k == 4
            iN1 = iNcon4; iN2 = iNcon6;
        elseif k == 5
            iN1 = iNcon5; iN2 = iNcon6;
        end
        iEL = iEL+1;
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
for i = 1:(lengthlevelZaxis-1) % ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬������б�����
        iNcon1 = iNO+1+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor_all+1;
        iNcon4 = iNcon1+lengthXYcoor_all+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor_all-1;
        iNcon4 = iNcon1+lengthXYcoor_all+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+stairColu_num*2-1;
    end
    if i < levelPstart1 % ���ǵײ���Ļǽ
        k_end = 1;
    else
        k_end = 2;
    end
    for k = 1:k_end % ����� % б��+ƽ��
        if k == 1
            iN1 = iNcon1; iN2 = iNcon3; iN3 = iNcon4; iN4 = iNcon2;
        elseif k == 2
            iN1 = iNcon3; iN2 = iNcon5; iN3 = iNcon6; iN4 = iNcon4;
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2, iN3, iN4,...    % �嵥Ԫ���ĸ��ڵ��
            ELE_iSUB, ELE_iWID);
    end
end
fprintf(fileID,'\n');

%% FLOORLOAD
fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');

LTNAME = ROOF; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0; % ¥�ݺ��� Ŀǰ�������������4/3.5�����0.��ÿ2.1һ��壬�ʿ��ܺ�4.2һ��7/3.5��ࡣ(������)
DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';

iNO = iNO_init; % ��ʼ��iNO
for i = 1:(lengthlevelZaxis-1) % ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬������б�����
        iNcon1 = iNO+1+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor_all+1;
        iNcon4 = iNcon1+lengthXYcoor_all+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor_all*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor_all-1;
        iNcon4 = iNcon1+lengthXYcoor_all+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+stairColu_num*2-1;
    end
    if i < levelPstart1 % ���ǵײ���Ļǽ
        k_end = 1;
    else
        k_end = 2;
    end
    for k = 1:k_end % ����� % б��+ƽ��
        if k == 1
            iN1 = iNcon1; iN2 = iNcon3; iN3 = iNcon4; iN4 = iNcon2;
        elseif k == 2
            iN1 = iNcon3; iN2 = iNcon5; iN3 = iNcon6; iN4 = iNcon4;
        end
        fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
            LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
            iN1, iN2, iN3, iN4);
    end
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;
end