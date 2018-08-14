%% function
% MGT stairs
%
% Xu Yi, 25th March 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_stair(fileID, iNO, iEL, CoC_stair, Deg_stair, facade_stair_R, levelZaxis_f, levelPstart1, stairColu_num, stairL, stairW, ~, ~, ROOF, Arc_itvl)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis_f = length(levelZaxis_f(:));

XYcoor_i = zeros(stairColu_num,2);   % ��ͲXoY�����1(X)��2(Y)�С�
XYcoor_o3 = zeros(lengthlevelZaxis_f,stairColu_num*3,2);	% ����*3����ͲXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

% ��Ͳ��
stairXYtemp1 = [stairL/2, stairW/2; -stairL/2, stairW/2]; % ԭʼ���꣬δ��ת��δת����������ϵ
stairXYtemp1Deg = coorDeg([0,0], [1,0], stairXYtemp1(1,:)); % ��Ͳ��һ����X��ĽǶ�
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
% Ļǽ��Ͳ
for j = levelPstart1:lengthlevelZaxis_f
    XYcoor_o_x = sqrt(facade_stair_R(j)^2 - (stairW/2)^2);  % Ļǽ��Ȧ�㵱y=��stairW/2ʱ��x�ľ���ֵ
    XYcoor_o_y = sqrt(facade_stair_R(j)^2 - (stairL/2)^2);  % Ļǽ��Ȧ�㵱x=��stairL/2ʱ��y�ľ���ֵ
    
    stairXYtemp = coorTrans([facade_stair_R(j), 0], -stairXYtemp1Deg); % б���˵����ꡣ
    
    stairXY_o_temp1 = [XYcoor_o_x, stairW/2; stairXYtemp; stairL/2, XYcoor_o_y];  % 3�����ģ��
    stairXY_o_temp2 = zeros(length(stairXY_o_temp1),2);
    for i = 1:length(stairXY_o_temp1)
        stairXY_o_temp2(i,:) = stairXY_o_temp1((length(stairXY_o_temp1)-i+1),:); % ��ʱ����
    end
    for i = 1:length(stairXY_o_temp1)
        stairXY_o_temp2(i,1) = -stairXY_o_temp1((length(stairXY_o_temp1)-i+1),1); % ��Y��Գ�
    end
    stairXY_o_temp12 = [stairXY_o_temp1; stairXY_o_temp2];
    stairXY_o_temp4 = zeros(length(stairXY_o_temp12),2);
    for i = 1:length(stairXY_o_temp12)
        stairXY_o_temp4(i,:) = stairXY_o_temp12((length(stairXY_o_temp12)-i+1),:); % ��ʱ����
    end
    for i = 1:length(stairXY_o_temp12)
        stairXY_o_temp4(i,2) = -stairXY_o_temp12((length(stairXY_o_temp12)-i+1),2); % ��X��Գ�
    end
    stairXY_o = [stairXY_o_temp1; stairXY_o_temp2; stairXY_o_temp4]; % Ļǽ��Ͳ��������
    
    for i = 1:length(stairXY_o)   % ����������
        [XYcoor_o3(j,i,:)] = coorTrans(stairXY_o(i,:), Deg_stair); % Ļǽ��Ͳ������ % ����ת����Ƕ�
    end
end

% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_stair(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_stair(2);
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_stair(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_stair(2);

lengthXYcoor_i = length(XYcoor_i); % ��Ͳÿ��ڵ���
lengthXYcoor_f = length(stairXY_o); % Ļǽÿ��ڵ���
lengthXYcoor_all = lengthXYcoor_i + lengthXYcoor_f;  % ÿ��ڵ�������

for i = 1:lengthlevelZaxis_f  % length(A(:)) A����Ԫ�ظ���
    % ��Ͳ����Ͳ % �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻��ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
    for j = 1:lengthXYcoor_i % �ڲ�4������
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis_f(i));
    end
    for j = 1:lengthXYcoor_f % �ⲿ12������
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis_f(i));
    end
end
% ��ֱ��������
iNO_main_end = iNO; % ���ڵ��յ㱸�ݣ�����ֱ�����ڵ���㱸�ݡ�
XY_Deg_num = zeros(lengthlevelZaxis_f,lengthXYcoor_f); % ��ֱ�����ĸ����ߵķָ��ڵ���
P_start = zeros(1,2); P_end = zeros(1,2);
for i = levelPstart1:lengthlevelZaxis_f
    for j = 1:lengthXYcoor_f % Ļǽ
        P_start(:) = XYcoor_o3(i,j,:);
        if j == lengthXYcoor_f
            P_end(:) = XYcoor_o3(i,1,:);
        else
            P_end(:) = XYcoor_o3(i,j+1,:);
        end
        [iNO, Deg_num] = MGT_arc_FE(fileID, iNO, levelZaxis_f(i), CoC_stair, P_start, P_end, Arc_itvl);
        XY_Deg_num(i,j) = Deg_num;
    end
end
% ��ֱ��������
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
for i = 1:(lengthlevelZaxis_f-1)	% length(A(:)) A����Ԫ�ظ���
    for j = 1:lengthXYcoor_i	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+j+lengthXYcoor_all*(i-1);
        iN2 = iN1+lengthXYcoor_all;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end

% ��Ͳ����iPRO = 1 ������1��
fprintf(fileID,'; Ļǽ��\n');
ELE_iPRO = 1;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart1:(lengthlevelZaxis_f-1)	% length(A(:)) A����Ԫ�ظ��� % levelPstart �ڼ��㿪ʼͣ�������¼��㿪��
    for j = 1:lengthXYcoor_f	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+(lengthXYcoor_i+j)+lengthXYcoor_all*(i-1); % ��������Ͳ��ͬ������ +stairN_num/2
        iN2 = iN1+lengthXYcoor_all;
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
fprintf(fileID,'; ¥������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = 1:lengthlevelZaxis_f	%
    for j = 1:lengthXYcoor_i	% ÿ����Ͳ�Ľڵ��� ����2~3(����Ϣƽ̨ʱ����platform.m)��4~1�ڵ㲻��(��Ϣƽ̨����Ҫ��)
        if j == 2 || j == 4
            % ���� 2~3��4~1�ڵ�
        else
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcoor_all*(i-1);
            iN2 = iN1+1;    %
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % ����Ԫ�������ڵ��
                ELE_ANGLE, ELE_iSUB);
        end
    end
end

fprintf(fileID,'; ¥��������\n');
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart1:lengthlevelZaxis_f % ¥�ݶ��Ǵ���Ͳ���
    for j = 1:lengthXYcoor_i % ÿ����Ͳ�Ľڵ���
        for k = 1:(lengthXYcoor_f/lengthXYcoor_i) % ÿ����Ͳ��Ӧ����Ͳ�ڵ���
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcoor_all*(i-1); % ����Ϊ��λ������¥�Ľڵ�(��Ͳ)
            iN2 = iNO+lengthXYcoor_i+(lengthXYcoor_f/lengthXYcoor_i)*(j-1)+k+lengthXYcoor_all*(i-1);    % �鵽Ļǽ��Ͳ��0����ٶ�λ�������
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
iNO_arc = iNO_main_end; % ��ʼ��
for i = levelPstart1:lengthlevelZaxis_f	% ����������Ԫ��ͬ������ԪΪi-1;
    for j = 1:lengthXYcoor_f	% ÿ����Ͳ�Ľڵ���
        iN1_bkp = iNO+lengthXYcoor_i+j+lengthXYcoor_all*(i-1); %
        if j ~= lengthXYcoor_f
            iN2_bkp = iN1_bkp+1;
        else % j = lengthXYcoor_f ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2_bkp = iN1_bkp+1-lengthXYcoor_f;
        end
        
        if XY_Deg_num(i,j) == 1 % ���˴�δ������ֱ�����ָ�
            iEL = iEL+1;
            iN1 = iN1_bkp;
            iN2 = iN2_bkp;
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % ����Ԫ�������ڵ��
                ELE_ANGLE, ELE_iSUB);
        else
            for k = 1:XY_Deg_num(i,j)
                iEL = iEL+1;
                if k == 1
                    iNO_arc = iNO_arc+1;
                    iN1 = iN1_bkp;
                    iN2 = iNO_arc;
                elseif k == XY_Deg_num(i,j)
                    iN1 = iNO_arc;
                    iN2 = iN2_bkp;
                else
                    iN1 = iNO_arc;
                    iNO_arc = iNO_arc+1;
                    iN2 = iNO_arc;
                end
                fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                    iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                    iN1, iN2,...    % ����Ԫ�������ڵ��
                    ELE_ANGLE, ELE_iSUB);
            end
        end
    end
end
fprintf(fileID,'\n');

%%
iEL_end = iEL;

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % ��ʼ��iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+lengthXYcoor_i);
CONSTRAINT = 111111; % 6�����ɶ�ȫԼ��
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end