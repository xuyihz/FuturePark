%% function
% MGT elevator
%
% Xu Yi, 24th April 2018
% Xu Yi, 24th April 2018, revised

%%
function [iNO_end, iEL_end] = MGT_elevator(fileID, iNO, iEL, CoC_elevator, Deg_elevator, facade_ele_R, levelZaxis_f, levelPstart1, elevatorColu_num, elevatorR, ~, ~, ROOF, Arc_itvl)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis_f = length(levelZaxis_f(:));

XYcoor_i = zeros(elevatorColu_num,2);   % 8����(�������м䲻������1�����Բ�ĵ�)��ͲXoY�����1(X)��2(Y)�С�
elevatorColu_o_num = elevatorColu_num+4; % ��Ͳ8����,��Ͳ12����
XYcoor_o3 = zeros(lengthlevelZaxis_f,elevatorColu_o_num,2);	% 12���㣬��ͲXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

% ��Ͳ��
elevatorXYtemp = elevatorR/sqrt(2); % ������
elevatorXY = [elevatorXYtemp,elevatorXYtemp; -elevatorXYtemp,elevatorXYtemp;...
    -elevatorXYtemp,-elevatorXYtemp; elevatorXYtemp,-elevatorXYtemp;...
    elevatorXYtemp,350; -elevatorXYtemp,350; 0,elevatorXYtemp; 0,350]; % ԭʼ���꣬δ��ת��δת����������ϵ
for i = 1:elevatorColu_num   % ����������
    XYcoor_i(i,:) = coorTrans(elevatorXY(i,:), Deg_elevator); % ��Ͳ
end
% ��ͲelevatorXY2 % ��Ͳ��Ҫ����Ļǽ���Ƥ���߶�λ % �����ظ߶ȵ�ѭ��(��������)
for j = levelPstart1:lengthlevelZaxis_f
    elevatorXY_o_temp = sqrt(facade_ele_R(j)^2 - elevatorXYtemp^2); %
    elevatorXY_o_Rtemp = sqrt(facade_ele_R(j)^2/2); % 2/5/8/11 �������Ľ�
    elevatorXY_o = [elevatorXY_o_temp,elevatorXYtemp; elevatorXY_o_Rtemp,elevatorXY_o_Rtemp; elevatorXYtemp,elevatorXY_o_temp;...
        -elevatorXYtemp,elevatorXY_o_temp; -elevatorXY_o_Rtemp,elevatorXY_o_Rtemp; -elevatorXY_o_temp,elevatorXYtemp;...
        -elevatorXY_o_temp,-elevatorXYtemp; -elevatorXY_o_Rtemp,-elevatorXY_o_Rtemp; -elevatorXYtemp,-elevatorXY_o_temp;...
        elevatorXYtemp,-elevatorXY_o_temp; elevatorXY_o_Rtemp,-elevatorXY_o_Rtemp; elevatorXY_o_temp,-elevatorXYtemp];
    
    for i = 1:elevatorColu_o_num   % ����������
        [XYcoor_o3(j,i,:)] = coorTrans(elevatorXY_o(i,:), Deg_elevator); % ��Ͳ
    end
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_i(:,1) = XYcoor_i(:,1) + CoC_elevator(1);
XYcoor_i(:,2) = XYcoor_i(:,2) + CoC_elevator(2);
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_elevator(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_elevator(2);

lengthXYcoor_i = length(XYcoor_i); % ��Ͳÿ��ڵ���
lengthXYcoor_f = length(elevatorXY_o); % Ļǽÿ��ڵ���
lengthXYcoor_all = lengthXYcoor_i + lengthXYcoor_f;  % ÿ��ڵ�������

for i = 1:lengthlevelZaxis_f  % length(A(:)) A����Ԫ�ظ���
    % ��Ͳ����Ͳ % �ڵ��Ź�����ÿ����Ͳ����ÿ����Ͳ�����µ��ϡ�
    for j = 1:elevatorColu_num % �ڲ�7������(��һ��������)
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor_i(j,1),XYcoor_i(j,2),levelZaxis_f(i));
    end
    for j = 1:elevatorColu_o_num % �ⲿ12������
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
        [iNO, Deg_num] = MGT_arc_FE(fileID, iNO, levelZaxis_f(i), CoC_elevator, P_start, P_end, Arc_itvl);
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
        if j==7 || j==8 % �����м�ڵ㲻����
        else
            iEL = iEL+1;
            iN1 = iNO+j+lengthXYcoor_all*(i-1);
            iN2 = iN1+lengthXYcoor_all;
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
for i = levelPstart1:(lengthlevelZaxis_f-1)	% length(A(:)) A����Ԫ�ظ��� % levelPstart �ڼ��㿪ʼͣ�������¼��㿪��
    for j = 1:lengthXYcoor_f	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO+(lengthXYcoor_i+j)+lengthXYcoor_all*(i-1); % ��������Ͳ��ͬ������ +elevatorColu_num
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
fprintf(fileID,'; ����Ͳ��Ͳ����\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
iN_table = [1,2; 2,3; 3,4; 4,1; 5,6; 7,8]; % ������ % ÿ��Ϊ�����˽ڵ��
for i = 1:lengthlevelZaxis_f	% ����i��ʼΪ2.�����㿪ʼ�С�
    for k = 1:length(iN_table) % ����6��
        iEL = iEL+1;
        iN1 = iNO + iN_table(k,1) + lengthXYcoor_all*(i-1);
        iN2 = iNO + iN_table(k,2) + lengthXYcoor_all*(i-1);
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % ����Ԫ�������ڵ��
            ELE_ANGLE, ELE_iSUB);
    end
end

% ��������iPRO = 3 ������3��
fprintf(fileID,'; ������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
for i = levelPstart1:lengthlevelZaxis_f	%
    for j = 1:lengthXYcoor_i/2 % ��Ͳ4����
        for k = 1:3 % ��Ͳÿ��������������
            iEL = iEL+1;
            iN1 = iNO + j + lengthXYcoor_all*(i-1); % ��Ͳ4���ǵ�
            iN2 = iNO + k + 3*(j-1) + lengthXYcoor_i + lengthXYcoor_all*(i-1); % ��Ͳ8����
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % ����Ԫ�������ڵ��
                ELE_ANGLE, ELE_iSUB);
        end
    end
end

% ��������iPRO = 4 ������4��
fprintf(fileID,'; ����Ͳ�⻷����\n');
ELE_iPRO = 4;
iNO = iNO_init; % ��ʼ��iNO
% Ļǽ�⻷��
iNO_arc = iNO_main_end; % ��ʼ��
for i = levelPstart1:lengthlevelZaxis_f	%
    for j = 1:lengthXYcoor_f
        iN1_bkp = iNO+j+lengthXYcoor_i+lengthXYcoor_all*(i-1); %
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
% 1~6 ���7���ڵ�û�������ӣ�����������7
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+lengthXYcoor_i-2);
CONSTRAINT = 111111; % 6�����ɶ�ȫԼ��
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end