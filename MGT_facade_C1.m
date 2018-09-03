%% function
% MGT tower 1 facade
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_facade_C1(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, facade_tower_R, levelZaxis, levelPstart1, iNO_towerC_init, Arc_itvl)
%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

iNO_init = iNO;
lengthlevelZaxis = length(levelZaxis(:));
XYcoor_o3 = zeros(lengthlevelZaxis,car_num*2,2);	% ��ĻǽXoY�����1(X)��2(Y)�С�ע����������ά���顣(Z����Ļǽ�б仯)

car_num2pi = 2*pi/car_num;  % speed up

XYcoor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % ����Y��ģ����Ͳһ�� X
XYcoor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y

for j = levelPstart1:lengthlevelZaxis
    XYcoor_o_1(1,1) = sqrt(facade_tower_R(j)^2 - XYcoor_i_1(1,2)^2);        % ����Y��ģ����Ͳһ�� X1 ע����Ͳ16���㲢���ǵȽǶȵȷ֡�
    XYcoor_o_1(1,2) = XYcoor_i_1(1,2);                                % Y1
    XYcoor_o_1(2,:) = coorMir(XYcoor_o_1(1,:), [0,0], XYcoor_i_1);     % X2,Y2
    for i = 0:(car_num-1)   % ���������� % ��ת�ֲ��Ƕ�+����Ƕ�
        [XYcoor_o3(j,i*2+1,:)] = coorTrans(XYcoor_o_1(1,:), -car_num2pi*i+Deg_tower); % ��Ͳ������1
        [XYcoor_o3(j,i*2+2,:)] = coorTrans(XYcoor_o_1(2,:), -car_num2pi*i+Deg_tower); % ��Ͳ������2
    end
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_o3(:,:,1) = XYcoor_o3(:,:,1) + CoC_tower(1);
XYcoor_o3(:,:,2) = XYcoor_o3(:,:,2) + CoC_tower(2);

lengthXYcoor_f = car_num*2;  % ÿ��ڵ���

for i = 1:lengthlevelZaxis  % length(A(:)) A����Ԫ�ظ���
    for j = 1:lengthXYcoor_f % Ļǽ
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_o3(i,j,1),XYcoor_o3(i,j,2),levelZaxis(i));   % ��Ͳ X & Y
    end
end
% ��ֱ��������
iNO_main_end = iNO; % ���ڵ��յ㱸�ݣ�����ֱ�����ڵ���㱸�ݡ�
XY_Deg_num = zeros(lengthlevelZaxis,lengthXYcoor_f); % ��ֱ�����ĸ����ߵķָ��ڵ���
P_start = zeros(1,2); P_end = zeros(1,2); % ÿ�����ߵ������յ�
for i = levelPstart1:lengthlevelZaxis
    for j = 1:lengthXYcoor_f % Ļǽ
        P_start(:) = XYcoor_o3(i,j,:);
        if j == lengthXYcoor_f
            P_end(:) = XYcoor_o3(i,1,:);
        else
            P_end(:) = XYcoor_o3(i,j+1,:);
        end
        [iNO, Deg_num] = MGT_arc_FE(fileID, iNO, levelZaxis(i), CoC_tower, P_start, P_end, Arc_itvl);
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

% ����������iPRO = 3 ������3��
fprintf(fileID,'; ��������\n');
ELE_iPRO = 3;
iNO = iNO_init; % ��ʼ��iNO
iNO_towerC = iNO_towerC_init; % ��ʼ��iNO_towerC
for i = levelPstart1:lengthlevelZaxis	% ������Ͳ���
    for j = 1:car_num*2	% ÿ����Ͳ�Ľڵ���
        iEL = iEL+1;
        iN1 = iNO_towerC+car_num+j+(car_num+car_num*2)*(i-1); % ����Ϊ��λ������¥�Ľڵ�(��Ͳ)
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
iNO_arc = iNO_main_end; % ��ʼ��
for i = levelPstart1:lengthlevelZaxis	% ����������Ԫ��ͬ������ԪΪi-1;
    for j = 1:lengthXYcoor_f	% ÿ����Ͳ�Ľڵ���
        iN1_bkp = iNO+j+lengthXYcoor_f*(i-1); % ÿ�����ߵ����
        if j ~= lengthXYcoor_f
            iN2_bkp = iN1_bkp+1; % ÿ�����ߵ��յ�
        else % j = lengthXYcoor_f ʱ�� ���ӵ��Ǳ����ĵ�һ���㣬�������ϲ��ڻ��ĵ�һ���㡣
            iN2_bkp = iN1_bkp+1-lengthXYcoor_f;
        end
        
        if XY_Deg_num(i,j) == 1 % ���˴�δ������ֱ�����ָ�������������С��Arc_itvl
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
                if k == 1 % ���ߵĵ�һС��
                    iNO_arc = iNO_arc+1;
                    iN1 = iN1_bkp;
                    iN2 = iNO_arc;
                elseif k == XY_Deg_num(i,j) % ���ߵ����һС��
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
iEL_end = iEL;
fprintf(fileID,'\n');
end