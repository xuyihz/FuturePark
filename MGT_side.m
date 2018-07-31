%% function
% MGT side tower
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_side(fileID, iNO, iEL, CoC_side, facade_side_R, levelZaxis, levelPstart, Roof_boundary, ~, ~, ROOF, tower_num)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;

switch tower_num % �����
    case 7
        sideRadius = 3300;
        Corner_coor = Roof_boundary(4,:);
        
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, sideRadius, Corner_coor, Roof_boundary(3,:)); %��ߵ�һ�� % [X, Y, Len] = coorLxC_sym(C0, R, P1, P2);
        if X_temp(1) < X_temp(2)
            XYcoor_1 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_1 = [X_temp(2),Y_temp(2)];
        end
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, sideRadius, Corner_coor, Roof_boundary(5,:)); %�ұߵ�һ�� % [X, Y, Len] = coorLxC_sym(C0, R, P1, P2);
        if X_temp(1) > X_temp(2)
            XYcoor_4 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_4 = [X_temp(2),Y_temp(2)];
        end
        Deg = coorDeg(Corner_coor, XYcoor_1, XYcoor_4);
        Temp_coor_trans = coorTransLoc(Corner_coor, XYcoor_1, -Deg/3);
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, sideRadius, Corner_coor, Temp_coor_trans);
        if X_temp(1) < X_temp(2)
            XYcoor_2 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_2 = [X_temp(2),Y_temp(2)];
        end
        Temp_coor_trans = coorTransLoc(Corner_coor, XYcoor_1, -Deg/3*2);
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, sideRadius, Corner_coor, Temp_coor_trans);
        if X_temp(1) < X_temp(2)
            XYcoor_3 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_3 = [X_temp(2),Y_temp(2)];
        end
        
        XYcoor_i = [XYcoor_1; XYcoor_2; XYcoor_3; XYcoor_4]; % ��ͲXoY�����1(X)��2(Y)�С� % ���ǵ�
        sideColu_num = length(XYcoor_i);
        
        lengthlevelZaxis = length(levelZaxis(:));
        XYcoor_o3 = zeros(lengthlevelZaxis,sideColu_num,2); % ���ƤXoY�����1(X)��2(Y)�С�
        for i = 1:(levelPstart-1)
            for j = 1:sideColu_num
                XYcoor_o3(i,j,:) = CoC_side;
            end
        end
        for i = levelPstart:lengthlevelZaxis % levelPstart���¶�Ϊ0
            for j = 1:sideColu_num
                [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, facade_side_R(i), Corner_coor, XYcoor_i(j,:));
                if Y_temp(1) < Y_temp(2)
                    XYcoor_temp = [X_temp(1),Y_temp(1)];
                else
                    XYcoor_temp = [X_temp(2),Y_temp(2)];
                end
                XYcoor_o3(i,j,:) = XYcoor_temp;
            end
        end
    case 10
        sideRadius = 2500;
        Corner_coor = Roof_boundary(2,:);
        
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, sideRadius, Corner_coor, Roof_boundary(1,:)); %��ߵ�һ�� % [X, Y, Len] = coorLxC_sym(C0, R, P1, P2);
        if X_temp(1) > X_temp(2)
            XYcoor_1 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_1 = [X_temp(2),Y_temp(2)];
        end
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, sideRadius, Corner_coor, Roof_boundary(3,:)); %�ұߵ�һ�� % [X, Y, Len] = coorLxC_sym(C0, R, P1, P2);
        if X_temp(1) > X_temp(2)
            XYcoor_3 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_3 = [X_temp(2),Y_temp(2)];
        end
        Deg = coorDeg(Corner_coor, XYcoor_1, XYcoor_3);
        Temp_coor_trans = coorTransLoc(Corner_coor, XYcoor_1, -Deg/2);
        [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, sideRadius, Corner_coor, Temp_coor_trans);
        if X_temp(1) > X_temp(2)
            XYcoor_2 = [X_temp(1),Y_temp(1)];
        else
            XYcoor_2 = [X_temp(2),Y_temp(2)];
        end
        
        XYcoor_i = [XYcoor_1; XYcoor_2; XYcoor_3]; % ��ͲXoY�����1(X)��2(Y)�С� % ���ǵ�
        sideColu_num = length(XYcoor_i);
        
        lengthlevelZaxis = length(levelZaxis(:));
        XYcoor_o3 = zeros(lengthlevelZaxis,sideColu_num,2); % ���ƤXoY�����1(X)��2(Y)�С�
        for i = 1:(levelPstart-1)
            for j = 1:sideColu_num
                XYcoor_o3(i,j,:) = CoC_side;
            end
        end
        for i = levelPstart:lengthlevelZaxis % levelPstart���¶�Ϊ0
            for j = 1:sideColu_num
                [X_temp, Y_temp, ~] = coorLxC_sym(CoC_side, facade_side_R(i), Corner_coor, XYcoor_i(j,:));
                if X_temp(1) > X_temp(2)
                    XYcoor_temp = [X_temp(1),Y_temp(1)];
                else
                    XYcoor_temp = [X_temp(2),Y_temp(2)];
                end
                XYcoor_o3(i,j,:) = XYcoor_temp;
            end
        end
end

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for k = 1:3 % �ǵ㣬��Ͳ����Ͳ
        if k == 1
            iNO = iNO+1;
            fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                iNO,Corner_coor(1),Corner_coor(2),levelZaxis(i));
        elseif k == 2 % �ڵ��Ź�����ʱ�롣
            for j = 1:sideColu_num % �ڲ�����
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis(i));
            end
        else
            for j = 1:sideColu_num % �ⲿ����
                iNO = iNO+1;
                fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
                    iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));
            end
        end
    end
end
lengthXYcoor2 = 1 + sideColu_num*2;  % ÿ��Ľڵ��������нǵ�1���㣬�ڲ�4���㣬�ⲿ4���㡣
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
    for j = 1:(1+sideColu_num)	% ÿ����Ͳ�Ľڵ��� 5
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
    for j = 1:sideColu_num	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+(1+sideColu_num+j)+lengthXYcoor2*(i-1); % ��������Ͳ��ͬ������ 1+sideColu_num
        iN2 = iN1+lengthXYcoor2;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame)
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% ������iPRO = 3 ������3��
fprintf(fileID,'; ��������״����\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis	% 
    for j = 1:sideColu_num
        iEL = iEL+1;
        iN1 = iNO+1+lengthXYcoor2*(i-1); % �ǵ�
        iN2 = iN1+j;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'; ��������״������\n');
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart:lengthlevelZaxis	% 
    for j = 1:sideColu_num
        iEL = iEL+1;
        iN1 = iNO+1+j+lengthXYcoor2*(i-1); % �ڲ���
        iN2 = iN1+sideColu_num;
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
% �ڻ���
fprintf(fileID,';   �ڻ���\n');
for i = 1:lengthlevelZaxis	% 
    for j = 1:(sideColu_num-1)
        iEL = iEL+1;
        iN1 = iNO+1+j+lengthXYcoor2*(i-1); % �ڲ���
        iN2 = iN1+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
% �⻷��
fprintf(fileID,';   �⻷��\n');
for i = levelPstart:lengthlevelZaxis	% 
    for j = 1:(sideColu_num-1)
        iEL = iEL+1;
        iN1 = iNO+1+sideColu_num+j+lengthXYcoor2*(i-1); % �ڲ���
        iN2 = iN1+1;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');
iEL_end = iEL;

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % ��ʼ��iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+1+sideColu_num);
CONSTRAINT = 111111; % 6�����ɶ�ȫԼ��
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end