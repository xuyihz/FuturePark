%% function
% MGT side tower
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_side(fileID, iNO, iEL, CoC_side, levelZaxis, levelPstart, Roof_boundary, ~, ~, ROOF, tower_num)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;

switch tower_num % �����
    case 7
        facade_stair_R = [zeros(6,1); 6820; 5839; 5162; 4734; 4529; 4534; 4750; 5189; 5880; 6878; 8473; 11040; 16200];
        Corner_coor = Roof_boundary(4,:);
        
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, 3250, Corner_coor, Roof_boundary(3,:)); %��ߵ�һ�� % [X, Y, Len] = coorLxC_sym(C0, R, P1, P2);
        if X_temp(1) < X_temp(2)
            XYcoor_1 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_1 = [X_temp(2),Y_temp(2)];
        end
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, 3250, Corner_coor, Roof_boundary(5,:)); %��ߵ�һ�� % [X, Y, Len] = coorLxC_sym(C0, R, P1, P2);
        if X_temp(1) > X_temp(2)
            XYcoor_4 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_4 = [X_temp(2),Y_temp(2)];
        end
        Deg = coorDeg(Corner_coor, XYcoor_1, XYcoor_4);
        Temp_coor = XYcoor_1 - Corner_coor;
        Temp_coor_trans = coorTrans(Temp_coor, -Deg/3);
        
    case 10
        facade_stair_R = [zeros(6,1); 7594; 6144; 5151; 4495; 4120; 4000; 4128; 4511; 5176; 6182; 7830; 10527; 16000];
        Corner_coor = Roof_boundary(2,:);
end

XYcoor_i = zeros(sideColu_num,2);   % ��ͲXoY�����1(X)��2(Y)�С� % ���ǵ�
XYcoor_o = zeros(sideColu_num,2);	% ���ƤXoY�����1(X)��2(Y)�С�

%
% �ֲ�����ϵ ת���� ��������ϵ

lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for k = 1:3 % �ǵ㣬��Ͳ����Ͳ
        if k == 1
            iNO = iNO+1;
            fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                iNO,Corner_coor(1),Corner_coor(2),levelZaxis(i));
        elseif k == 1 % �ڵ��Ź�����ʱ�롣
            for j = 1:sideColu_num % �ڲ�����
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
            end
        else
            for j = 1:sideColu_num % �ⲿ����
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_o(j,1),XYcoor_o(j,2),levelZaxis(i));
            end
        end
    end
end
lengthXYcoor2 = sideColu_num+sideColu_num*2;  % ÿ��Ľڵ����������ڲ�4���㣬�ⲿ8���㡣
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
    for j = 1:sideColu_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor2*(i-1);
        iN2 = iN1+lengthXYcoor2;
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
    for j = 1:sideColu_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+(sideColu_num+j)+lengthXYcoor2*(i-1); % ��������Ͳ��ͬ������ +stairN_num/2
        iN2 = iN1+lengthXYcoor2;
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

% ������iPRO = 3 ������3��
fprintf(fileID,'; ¥�ݳ�������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = 1:(lengthlevelZaxis-1)	% ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬������б�����
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor2+1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor2-1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+sideColu_num*2-1;
    end
    if i < levelPstart % ���ǵײ���Ļǽ
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
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+1;
    end
    if i < levelPstart % ���ǵײ���Ļǽ
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
for i = levelPstart:(lengthlevelZaxis-1)	% ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬������б�����
        iNcon1 = iNO+6+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+5;
        iNcon3 = iNcon1+lengthXYcoor2+1;
        iNcon4 = iNcon2+lengthXYcoor2-1;
        iNcon5 = iNcon3+1;
        iNcon6 = iNcon4-1;
    else % ż����
        iNcon1 = iNO+7+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor2-1;
        iNcon4 = iNcon2+lengthXYcoor2+1;
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
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor2+1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor2-1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+sideColu_num*2-1;
    end
    if i < levelPstart % ���ǵײ���Ļǽ
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
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
        iNcon3 = iNcon1+lengthXYcoor2+1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+6;
        iNcon6 = iNcon5+1;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+1;
        iNcon3 = iNcon1+lengthXYcoor2-1;
        iNcon4 = iNcon1+lengthXYcoor2+2;
        iNcon5 = iNcon3+4;
        iNcon6 = iNcon5+sideColu_num*2-1;
    end
    if i < levelPstart % ���ǵײ���Ļǽ
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

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % ��ʼ��iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+sideColu_num);
CONSTRAINT = 111111; % 6�����ɶ�ȫԼ��
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end