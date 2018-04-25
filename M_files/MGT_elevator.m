%% function
% MGT elevator
%
% Xu Yi, 24th April 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_elevator(fileID, iNO, iEL, CoC_elevator, Deg_elevator, levelZaxis, levelPstart, elevatorColu_num, ~, ~, ROOF)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;

elevatorColu_o_num = elevatorColu_num+2;
XYcor_i = zeros(elevatorColu_num,2);   % 7�� ��ͲXoY�����1(X)��2(Y)�С�
XYcor_o = zeros(elevatorColu_o_num,2); % 10�� ���ƤXoY�����1(X)��2(Y)�С�

ele_width = 2500; ele_shift = 400; ele_depth = ele_width + ele_shift; 
str_length = 2200; str_width = 2800;
elevatorXY = [-ele_width, ele_depth; ele_width, ele_depth; -ele_width, ele_shift; ele_width, ele_shift;...
                -str_length, -str_width; str_length, -str_width; 0, ele_depth; 0, ele_shift]; % ԭʼ���꣬δ��ת��δת����������ϵ
for i = 1:elevatorColu_num   % ����������
    [XYcor_i(i,1), XYcor_i(i,2)] = coorTrans(elevatorXY(i,1), elevatorXY(i,2), Deg_elevator); % ��Ͳ
end
% ��Ͳ����stairXY2 �ݰ��뾶6000��λ
ele2_depth = 5450; ele2_width1 = 5200; ele2_width2 = 6000; ele2_width3 = 5300;
ele2_strX = str_length*2 - ele_width; ele2_strY = str_width*2 + ele_shift;
elevatorXY2 = [-ele_width, ele2_depth; ele_width, ele2_depth; -ele2_width1, ele_depth; ele2_width1, ele_depth;...
                -ele2_width2, ele_shift; ele2_width2, ele_shift; -ele2_width3, -str_width; ele2_width3, -str_width;...
                -ele2_strX, -ele2_strY; ele2_strX, -ele2_strY];
for i = 1:elevatorColu_o_num   % ����������
    [XYcor_o(i,1), XYcor_o(i,2)] = coorTrans(elevatorXY2(i,1), elevatorXY2(i,2), Deg_elevator); % ��Ͳ
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcor_i(:,1) = XYcor_i(:,1) + CoC_elevator(1);
XYcor_i(:,2) = XYcor_i(:,2) + CoC_elevator(2);
XYcor_o(:,1) = XYcor_o(:,1) + CoC_elevator(1);
XYcor_o(:,2) = XYcor_o(:,2) + CoC_elevator(2);
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for k = 1:2 % ��Ͳ����Ͳ
        if k == 1 % �ڵ��Ź�����ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
            for j = 1:elevatorColu_num % �ڲ�7������
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcor_i(j,1),XYcor_i(j,2),levelZaxis(i));
            end
        else
            for j = 1:elevatorColu_o_num % �ⲿ10������
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcor_o(j,1),XYcor_o(j,2),levelZaxis(i));
            end
        end
    end
end
lengthXYcor2 = elevatorColu_num+elevatorColu_o_num;  % ÿ��Ľڵ����������ڲ�47���㣬�ⲿ10���㡣
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
    for j = 1:elevatorColu_num	% ÿ����Ͳ�Ľڵ���
        if j == 7   % �����м��Ǹ��ڵ㲻����
        else
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcor2*(i-1);
            iN2 = iN1+lengthXYcor2;
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % ����Ԫ�������ڵ��
                ELE_ANGLE, ELE_iSUB);
        end
    end
end

% ��Ͳ����iPRO = 1 ������1��
fprintf(fileID,'; ��Ͳ��\n');
ELE_iPRO = 1;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart:(lengthlevelZaxis-1)	% length(A(:)) A����Ԫ�ظ��� % levelPstart �ڼ��㿪ʼͣ�������¼��㿪��
    for j = 1:elevatorColu_o_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+(elevatorColu_num+j)+lengthXYcor2*(i-1); % ��������Ͳ��ͬ������ +elevatorColu_num
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

% ������iPRO = 3 ������3��
fprintf(fileID,'; ��������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = 2:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1; ����i��ʼΪ2.�����㿪ʼ�С�
    iNcon = iNO+lengthXYcor2*(i-1);
    for k = 1:7 % ����7��
        switch k
            case 1
                iN1 = iNcon+1; iN2 = iNcon+3;
            case 2
                iN1 = iNcon+7; iN2 = iNcon+8;
            case 3
                iN1 = iNcon+2; iN2 = iNcon+4;
            case 4
                iN1 = iNcon+1; iN2 = iNcon+7;
            case 5
                iN1 = iNcon+7; iN2 = iNcon+2;
            case 6
                iN1 = iNcon+3; iN2 = iNcon+8;
            case 7
                iN1 = iNcon+8; iN2 = iNcon+4;
        end        
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end

% ������iPRO = 3 ������3��
fprintf(fileID,'; ¥�ݳ�������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = 1:(lengthlevelZaxis-1)	% ������б�Σ�������Ҫ-1
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬������б�����
        iN1 = iNO+5+lengthXYcor2*(i-1); % �ݶ��ڵ�5���
        iN2 = iN1+lengthXYcor2+1;
    else % ż����
        iN1 = iNO+6+lengthXYcor2*(i-1);
        iN2 = iN1+lengthXYcor2-1;
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % ����Ԫ�������ڵ��
        ELE_ANGLE, ELE_iSUB);
end
fprintf(fileID,'; ¥�ݿ�������\n');
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1 % ÿ��һ���ᴩ��������
    if rem(i,2) ~= 0    % ������ % ���Ƶ㣬��������Ͳ�������
        iN1 = iNO+3+lengthXYcor2*(i-1);
        iN2 = iN1+2;
    else % ż����
        iN1 = iNO+4+lengthXYcor2*(i-1);
        iN2 = iN1+2;
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % ����Ԫ�������ڵ��
        ELE_ANGLE, ELE_iSUB);
end

% ��������iPRO = 3 ������3��
fprintf(fileID,'; ������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart:lengthlevelZaxis	%
    iNcon = iNO+lengthXYcor2*(i-1);
    iNcon_o = iNcon + elevatorColu_num; % ��Ͳ�ڵ����
    for k = 1:10 % ��������10��
        switch k
            case 1
                iN1 = iNcon+1; iN2 = iNcon_o+1;
            case 2
                iN1 = iNcon+1; iN2 = iNcon_o+3;
            case 3
                iN1 = iNcon+2; iN2 = iNcon_o+2;
            case 4
                iN1 = iNcon+2; iN2 = iNcon_o+4;
            case 5
                iN1 = iNcon+3; iN2 = iNcon_o+5;
            case 6
                iN1 = iNcon+4; iN2 = iNcon_o+6;
            case 7
                iN1 = iNcon+5; iN2 = iNcon_o+7;
            case 8
                iN1 = iNcon+5; iN2 = iNcon_o+9;
            case 9
                iN1 = iNcon+6; iN2 = iNcon_o+8;
            case 10
                iN1 = iNcon+6; iN2 = iNcon_o+10;
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
for i = levelPstart:lengthlevelZaxis	%
    iNcon_o = iNO+lengthXYcor2*(i-1) + elevatorColu_num; % ��Ͳ�ڵ����
    for k = 1:10 % ����10��
        switch k
            case 1
                iN1 = iNcon_o+1; iN2 = iNcon_o+2;
            case 2
                iN1 = iNcon_o+1; iN2 = iNcon_o+3;
            case 3
                iN1 = iNcon_o+2; iN2 = iNcon_o+4;
            case 4
                iN1 = iNcon_o+3; iN2 = iNcon_o+5;
            case 5
                iN1 = iNcon_o+4; iN2 = iNcon_o+6;
            case 6
                iN1 = iNcon_o+5; iN2 = iNcon_o+7;
            case 7
                iN1 = iNcon_o+6; iN2 = iNcon_o+8;
            case 8
                iN1 = iNcon_o+7; iN2 = iNcon_o+9;
            case 9
                iN1 = iNcon_o+8; iN2 = iNcon_o+10;
            case 10
                iN1 = iNcon_o+9; iN2 = iNcon_o+10;
        end
        iEL = iEL+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(planner) floor ���ڸ�4���㲻��һ��ƽ�棬���޷�ʩ��¥�����
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
        iN1 = iNO+3+lengthXYcor2*(i-1);
        iN2 = iN1+2;
        iN3 = iN1+3+lengthXYcor2;
        iN4 = iN1+1+lengthXYcor2;
    else % ż����
        iN1 = iNO+6+lengthXYcor2*(i-1);
        iN2 = iN1-2;
        iN3 = iN1-3+lengthXYcor2;
        iN4 = iN1-1+lengthXYcor2;
    end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2, iN3, iN4,...    % �嵥Ԫ���ĸ��ڵ��
        ELE_iSUB, ELE_iWID);
end
fprintf(fileID,'\n');

%% FLOORLOAD ���ڸ�4���㲻��һ��ƽ�棬���޷�ʩ��¥�����
% fprintf(fileID,'*FLOORLOAD    ; Floor Loads\n');
% fprintf(fileID,'; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn  ; iDIST=1,2\n; LTNAME, iDIST, DIR, bPROJ, DESC, GROUP, NODE1, ..., NODEn                                        ; iDIST=3,4\n; [iDIST] 1=One Way, 2=Two Way, 3=Polygon-Centroid, 4=Polygon-Length\n');
% 
% LTNAME = ROOF; iDIST = 2; ANGLE = 0; iSBEAM = 0; SBANG = 0; SBUW = 0; % ¥�ݺ��� Ŀǰ�������������4/3.5�����0.��ÿ2.1һ��壬�ʿ��ܺ�4.2һ��7/3.5��ࡣ(������)
% DIR = 'GZ'; bPROJ = 'NO'; DESC = ''; bEX = 'NO'; bAL = 'NO'; GROUP = '';
% 
% iNO = iNO_init; % ��ʼ��iNO
% for i = 1:(lengthlevelZaxis-1) % ������б�Σ�������Ҫ-1
%     if rem(i,2) ~= 0    % ������ % ���Ƶ㣬������б�����
%         iN1 = iNO+3+lengthXYcor2*(i-1);
%         iN2 = iN1+2;
%         iN3 = iN1+3+lengthXYcor2;
%         iN4 = iN1+1+lengthXYcor2;
%     else % ż����
%         iN1 = iNO+6+lengthXYcor2*(i-1);
%         iN2 = iN1-2;
%         iN3 = iN1-3+lengthXYcor2;
%         iN4 = iN1-1+lengthXYcor2;
%     end
%     fprintf(fileID,'   %s, %d, %d, %d, %d, %d, %s, %s, %s, %s, %s, %s, %d, %d, %d, %d\n',...
%         LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP,...
%         iN1, iN2, iN3, iN4);
% end
% fprintf(fileID,'\n');

%%
iEL_end = iEL;

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % ��ʼ��iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+elevatorColu_num);
CONSTRAINT = 111111; % 6�����ɶ�ȫԼ��
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end