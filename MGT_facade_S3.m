%% function
% MGT tower 3 facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_facade_S3(fileID, iNO, iEL, column_num, CoC_tower, Deg_tower, towerS_column_coor, facade_tower_R, levelZaxis, levelPstart, iNO_towerS_init)    %ע�������levelPstart��1x2����
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
levelPstart1 = levelPstart(1);
lengthlevelZaxis = length(levelZaxis(:));
xtemp = towerS_column_coor(1); ytemp = towerS_column_coor(2);
% XYcoor_i = [ xtemp,ytemp; -xtemp,ytemp; -xtemp,-ytemp; xtemp,-ytemp ]; % С��������
XYcoor_o3 = zeros(lengthlevelZaxis,column_num*2,2);	% ��ͲXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

Rec_R = sqrt( xtemp^2 + ytemp^2 ); % 4�������γɵľ��ζԽ��ߵ�һ��
for j = levelPstart1:lengthlevelZaxis
    Rtemp = facade_tower_R(j);
    XY_o_2x = xtemp / Rec_R * Rtemp;	% Ļǽ��2�� x��
    XY_o_2y = ytemp / Rec_R * Rtemp;  % Ļǽ��2�� y��
    
    XYcoor_o3(j,:,:) = [Rtemp, 0; XY_o_2x, XY_o_2y; 0, Rtemp; -XY_o_2x, XY_o_2y;...
        -Rtemp, 0; -XY_o_2x, -XY_o_2y; 0, -Rtemp; XY_o_2x, -XY_o_2y;];
    if Deg_tower ~= 0
        for i = 1:column_num*2
            XYcoor_o3(j,i,:) = coorTrans(XYcoor_o3(j,i,:), Deg_tower);
        end
    end
end

% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_tower(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_tower(2);

lengthXYcoor_f = column_num*2;  % Ļǽÿ��ڵ���
lengthlevelZaxis = length(levelZaxis(:));

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for j = 1:lengthXYcoor_f % Ļǽ
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));   % ��Ͳ X & Y
    end
end
iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_colu = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% ��Ͳ����iPRO = 1 ������1��
fprintf(fileID,'; Ļǽ��(������Ͳ��)\n');
ELE_iPRO = 1;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart1:(lengthlevelZaxis-1)	% length(A(:)) A����Ԫ�ظ��� % levelPstart �ڼ��㿪ʼͣ�������¼��㿪��
    for j = 1:lengthXYcoor_f	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_f*(i-1); % ����Ļǽ���������Ĳ�ͬ����Ļǽ�ڵ�Ϊ������ţ���һ��ֻ��lengthXYcoor_f���ڵ�
        iN2 = iN1+lengthXYcoor_f;
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
iNO_towerS = iNO_towerS_init; % ��ʼ��iNO_towerS
for i = levelPstart1:lengthlevelZaxis	% ������Ͳ���
    for j = 1:column_num	% ÿ����Ͳ�Ľڵ���
        for k = 1:3 % һ����Ͳ������������Ͳ������������
            iEL = iEL+1;
            iN1 = iNO_towerS+j+column_num*(i-1); % ����Ϊ��λ������¥�Ľڵ�(��Ͳ)
            iN2 = iNO+lengthXYcoor_f*(i-1)+(j-1)*2+k;    % �鵽Ļǽ��Ͳ��0����ٶ�λ�������
            if j == 4 && k == 3
                iN2 = iNO+lengthXYcoor_f*(i-1)+(j-1)*2+k-lengthXYcoor_f;
            end
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
% �⻷��
fprintf(fileID,';   Ļǽ�⻷��\n');
for i = levelPstart1:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1;
    for j = 1:column_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_f*(i-1); %
        if j ~= column_num*2
            iN2 = iN1+1;
        else % j = car_num*2 ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2 = iN1+1-column_num*2;
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