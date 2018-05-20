%% function
% MGT stairs
%
% Xu Yi, 25th March 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_stair(fileID, iNO, iEL, CoC_stair, Deg_stair, levelZaxis, levelPstart1, stairColu_num, stairL, stairW, stairB, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;

XYcoor_i = zeros(stairColu_num,2);   % ��ͲXoY�����1(X)��2(Y)�С�
XYcoor_o = zeros(stairColu_num*2,2); % ���ƤXoY�����1(X)��2(Y)�С�

stairXY = [stairL/2, stairW/2; -stairL/2, stairW/2; -stairL/2, -stairW/2; stairL/2, -stairW/2]; % ԭʼ���꣬δ��ת��δת����������ϵ
for i = 1:stairColu_num   % ����������
    [XYcoor_i(i,:)] = coorTrans(stairXY(i,:), Deg_stair); % ��Ͳ
end
% ��Ͳ����stairXY2 �ݶ�����Ͳ����stairB.
stairXY2 = [stairL/2+stairB, stairW/2; stairL/2, stairW/2+stairB; -stairL/2, stairW/2+stairB; -stairL/2-stairB, stairW/2;...
    -stairL/2-stairB, -stairW/2; -stairL/2, -stairW/2-stairB; stairL/2, -stairW/2-stairB; stairL/2+stairB, -stairW/2];
for i = 1:stairColu_num*2   % ����������
    [XYcoor_o(i,:)] = coorTrans(stairXY2(i,:), Deg_stair); % ��Ͳ
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_stair(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_stair(2);
XYcoor_o(:,1) = XYcoor_o(:,1) + CoC_stair(1);
XYcoor_o(:,2) = XYcoor_o(:,2) + CoC_stair(2);
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for k = 1:2 % ��Ͳ����Ͳ
        if k == 1 % �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻��ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
            for j = 1:stairColu_num % �ڲ�4������
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
            end
        else
            for j = 1:stairColu_num*2 % �ⲿ8������
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_o(j,1),XYcoor_o(j,2),levelZaxis(i));
            end
        end
    end
end
lengthXYcoor2 = stairColu_num+stairColu_num*2;  % ÿ��Ľڵ����������ڲ�4���㣬�ⲿ8���㡣
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
    for j = 1:stairColu_num	% ÿ����Ͳ�Ľڵ���
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
for i = levelPstart1:(lengthlevelZaxis-1)	% length(A(:)) A����Ԫ�ظ��� % levelPstart �ڼ��㿪ʼͣ�������¼��㿪��
    for j = 1:stairColu_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+(stairColu_num+j)+lengthXYcoor2*(i-1); % ��������Ͳ��ͬ������ +stairN_num/2
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
        iNcon1 = iNO+1+lengthXYcoor2*(i-1);
        iNcon2 = iNcon1+3;
    else % ż����
        iNcon1 = iNO+2+lengthXYcoor2*(i-1);
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

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % ��ʼ��iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+stairColu_num);
CONSTRAINT = 111111; % 6�����ɶ�ȫԼ��
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end