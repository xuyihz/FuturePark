%% function
% MGT stair facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_facade_stair(fileID, iNO, iEL, stairColu_num, CoC_stair, Deg_stair, stairL, stairW, levelZaxis, levelPstart1, iNO_stair_init, tower_num)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

% ����Ļǽÿ���������б仯������Ҫѭ��(��������)��������йأ�����Ƶ㶨λ�йأ����ݶ��̶�����ֵ(Ŀǰ��19��(����ƽ̨)����7����19��)�������������
switch tower_num % �����
    case {5, 6}
        facade_stair_R = [zeros(6,1); 6820; 5839; 5162; 4734; 4529; 4534; 4750; 5189; 5880; 6878; 8473; 11040; 16200];
    case {7, 9}
        facade_stair_R = [zeros(6,1); 7594; 6144; 5151; 4495; 4120; 4000; 4128; 4511; 5176; 6182; 7830; 10527; 16000];
    case 8
        facade_stair_R = [zeros(6,1); 9769; 8319; 7326; 6669; 6295; 6175; 6302; 6685; 7351; 8357; 10005; 12702; 18175];
    case 10
        facade_stair_R = [zeros(6,1); 6927; 5477; 4483; 3827; 3452; 3333; 3460; 3843; 4509; 5514; 7163; 9860; 15333];
end

iNO_init = iNO;
lengthlevelZaxis = length(levelZaxis(:));
XYcoor_i = zeros(stairColu_num,2);    % ��ͲXoY�����1(X)��2(Y)�С�
XYcoor_o3 = zeros(lengthlevelZaxis,stairColu_num*2,2);	% ��ͲXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

stairXY = [stairL/2, stairW/2; -stairL/2, stairW/2; -stairL/2, -stairW/2; stairL/2, -stairW/2]; % ԭʼ���꣬δ��ת��δת����������ϵ
for i = 1:stairColu_num   % ����������
    [XYcoor_i(i,1), XYcoor_i(i,2)] = coorTrans(stairXY(i,1), stairXY(i,2), Deg_stair); % ��Ͳ
end
% Ļǽ��Ͳ
for j = levelPstart1:lengthlevelZaxis
    XYcoor_o_x = sqrt(facade_stair_R(j)^2 - (stairW/2)^2);  % Ļǽ��Ȧ�㵱y=��stairW/2ʱ��x�ľ���ֵ
    XYcoor_o_y = sqrt(facade_stair_R(j)^2 - (stairL/2)^2);  % Ļǽ��Ȧ�㵱x=��stairL/2ʱ��y�ľ���ֵ
    stairXY2 = [XYcoor_o_x, stairW/2; stairL/2, XYcoor_o_y; -stairL/2, XYcoor_o_y; -XYcoor_o_x, stairW/2;...
        -XYcoor_o_x, -stairW/2; -stairL/2, -XYcoor_o_y; stairL/2, -XYcoor_o_y; XYcoor_o_x, -stairW/2];  % Ļǽ��Ͳ��������
    for i = 1:stairColu_num*2   % ����������
        [XYcoor_o3(j,i,1), XYcoor_o3(j,i,2)] = coorTrans(stairXY2(i,1), stairXY2(i,2), Deg_stair); % Ļǽ��Ͳ������ % ����ת����Ƕ�
    end
end

% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_stair(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_stair(2);
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_stair(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_stair(2);

lengthXYcoor2 = length(XYcoor_i(:))/2 + stairColu_num*2;  % ÿ��ڵ�������
lengthXYcoor_f = stairColu_num*2;  % Ļǽÿ��ڵ���

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for j = 1:(stairColu_num*2)
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));   % ��Ͳ X1 & Y1
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
    for j = 1:stairColu_num*2	% ÿ����Ͳ�Ľڵ���
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
iNO_stair = iNO_stair_init; % ��ʼ��iNO_stair
for i = levelPstart1:lengthlevelZaxis	% ¥�ݶ��Ǵ���Ͳ���
    for j = 1:stairColu_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO_stair+(stairColu_num+j)+lengthXYcoor2*(i-1); % ����Ϊ��λ������¥�Ľڵ�(��Ͳ)
        iN2 = iNO+lengthXYcoor_f*(i-1)+j;    % �鵽Ļǽ��Ͳ��0����ٶ�λ�������
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
% �⻷��
fprintf(fileID,';   Ļǽ�⻷��\n');
for i = levelPstart1:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1;
    for j = 1:stairColu_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_f*(i-1); % 
        if j ~= stairColu_num*2
            iN2 = iN1+1;
        else % j = car_num*2 ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2 = iN1+1-stairColu_num*2;
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