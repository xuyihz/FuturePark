%% function
% MGT towerC ramp Misc.
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_Misc(fileID, iNO, iEL, car_num, CoC_tower, Deg_tower, tube_innerR, ~, levelZaxis, levelPstart, Roof_boundary, ~,...
    CoC_towerS2, CoC_towerS3, CoC_elevator4, CoC_stair5, CoC_stair6, facade_tower2_R, facade_tower3_R, facade_ele4_R, facade_stair5_R, facade_stair6_R,...
    levelTaxis, levelTSaxis, levelSaxis_f)    % ע�������levelPstart��1x2����
%%
CoC = [CoC_tower; CoC_towerS2; CoC_towerS3; CoC_elevator4; CoC_stair5; CoC_stair6];                     % centre of tower1~6
facade_R = {[0,0]; facade_tower2_R; facade_tower3_R; facade_ele4_R; facade_stair5_R; facade_stair6_R};    % facade R of tower1~6
levelZ_f = {levelTaxis; levelTSaxis; levelTSaxis; levelSaxis_f; levelSaxis_f; levelSaxis_f};            % facade Z axis of tower1~6
f_boundary = {[0,0,0,0]; [Roof_boundary(1,:),Roof_boundary(2,:)]; [Roof_boundary(7,:),Roof_boundary(1,:)]; [Roof_boundary(7,:),Roof_boundary(1,:)];...
    [Roof_boundary(2,:),Roof_boundary(3,:),Roof_boundary(3,:),Roof_boundary(4,:)]; [Roof_boundary(4,:),Roof_boundary(5,:),Roof_boundary(5,:),Roof_boundary(6,:)]};  % 1~6��¥�Ӵ��ı���
%% NODE ������ҵ�����������µ�ͶӰ���Ľڵ�
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

% �����µ�ÿ���������б仯������Ҫѭ��(��������)��(��ʱδ���������仯) % ��17600��ʼ���¡����Ͻ���
ramp_inner_R = 12200; % �µ���Ȧ�뾶
ramp_outer_R = 14400; % �µ���Ȧ�뾶
aisle_width = 1800; % �ߵ����

iNO_init = iNO;
% levelPstart1 = levelPstart(1);
levelPstart2 = levelPstart(2);
% XYcoor_i = zeros(car_num,2);	% ��ͲXoY�����1(X)��2(Y)�С�
% XYcoor_o = zeros(car_num*2,2);	% ��ͲXoY�����1(X)��2(Y)�С�
ramp_point = 1;
XYcoor_ramp_i = zeros(ramp_point,car_num*2,2);	% �µ��ڲ�XoY�����1(X)��2(Y)�С�
XYcoor_ramp_o = zeros(ramp_point,car_num*2,2);	% �µ����XoY�����1(X)��2(Y)�С�
XYcoor_aisle = zeros(ramp_point,car_num*2,2);	% �ߵ����XoY�����1(X)��2(Y)�С�(�ߵ��ڲ༴Ϊ�µ����)
circle_num = 3; % ��������Ȧ

car_num2pi = 2*pi/car_num;  % speed up
Deg_tower = Deg_tower + pi/2; % �����µ���λ

XYcoor_i_1(1,1) = tube_innerR * cos(car_num2pi/2);   % ����Y��ģ����Ͳһ�� X % Y��ģ�鿴Parts�ļ����ڵ�Ͳ�ļ�
XYcoor_i_1(1,2) = tube_innerR * sin(car_num2pi/2);   % Y

XYcoor_ramp_i_1(1,1) = sqrt(ramp_inner_R^2 - XYcoor_i_1(1,2)^2);        % Y���µ��ڲ���һ�� X1 ע��16���㲢���ǵȽǶȵȷ֡�
XYcoor_ramp_i_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_ramp_i_1(2,:) = coorMir(XYcoor_ramp_i_1(1,:), [0,0], XYcoor_i_1);% ��Y�������߾��񣬵õ�Y��ģ�������֧��X2,Y2
XYcoor_ramp_o_1(1,1) = sqrt(ramp_outer_R^2 - XYcoor_i_1(1,2)^2);        % Y���µ������һ�� X1 ע��16���㲢���ǵȽǶȵȷ֡�
XYcoor_ramp_o_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_ramp_o_1(2,:) = coorMir(XYcoor_ramp_o_1(1,:), [0,0], XYcoor_i_1);% ��Y�������߾��񣬵õ�Y��ģ�������֧��X2,Y2
XYcoor_aisle_1(1,1) = sqrt((ramp_outer_R+aisle_width)^2 - XYcoor_i_1(1,2)^2);        % Y���ߵ������һ�� X1 ע��16���㲢���ǵȽǶȵȷ֡�
XYcoor_aisle_1(1,2) = XYcoor_i_1(1,2);                                 % Y1
XYcoor_aisle_1(2,:) = coorMir(XYcoor_aisle_1(1,:), [0,0], XYcoor_i_1);% ��Y�������߾��񣬵õ�Y��ģ�������֧��X2,Y2
for j = 1:ramp_point % �ݶ�1�������µ��������Ⱥ͸߶��йأ����ڻ�����ġ�
    for i = 0:(car_num-1)   % ���������� % ��ת�ֲ��Ƕ�+����Ƕ�
        [XYcoor_ramp_i(j,i*2+1,:)] = coorTrans(XYcoor_ramp_i_1(1,:), -car_num2pi*i+Deg_tower); % �µ��ڲ������1
        [XYcoor_ramp_i(j,i*2+2,:)] = coorTrans(XYcoor_ramp_i_1(2,:), -car_num2pi*i+Deg_tower); % �µ��ڲ������2
        [XYcoor_ramp_o(j,i*2+1,:)] = coorTrans(XYcoor_ramp_o_1(1,:), -car_num2pi*i+Deg_tower); % �µ���������1
        [XYcoor_ramp_o(j,i*2+2,:)] = coorTrans(XYcoor_ramp_o_1(2,:), -car_num2pi*i+Deg_tower); % �µ���������2
        [XYcoor_aisle(j,i*2+1,:)] = coorTrans(XYcoor_aisle_1(1,:), -car_num2pi*i+Deg_tower); % �µ���������1
        [XYcoor_aisle(j,i*2+2,:)] = coorTrans(XYcoor_aisle_1(2,:), -car_num2pi*i+Deg_tower); % �µ���������2
    end
end
% �ֲ�����ϵ ת���� ��������ϵ
XYcoor_ramp_i(:,:,1) = XYcoor_ramp_i(:,:,1) + CoC_tower(1);
XYcoor_ramp_i(:,:,2) = XYcoor_ramp_i(:,:,2) + CoC_tower(2);
XYcoor_ramp_o(:,:,1) = XYcoor_ramp_o(:,:,1) + CoC_tower(1);
XYcoor_ramp_o(:,:,2) = XYcoor_ramp_o(:,:,2) + CoC_tower(2);
XYcoor_aisle(:,:,1) = XYcoor_aisle(:,:,1) + CoC_tower(1);
XYcoor_aisle(:,:,2) = XYcoor_aisle(:,:,2) + CoC_tower(2);

[~,lengthXYcoor_ramp,~] = size(XYcoor_ramp_i);  % �µ�ÿ��ڵ���
for k = 1:2
    if k == 1
        levelZ = levelZaxis(levelPstart2);
    else
        levelZ = levelZaxis(end);
    end
    i = ramp_point;
    for j = 1:lengthXYcoor_ramp  % ��ҵ�㡢����� �µ�x������ͶӰ�ڵ�
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_ramp_i(i,j,1),XYcoor_ramp_i(i,j,2),levelZ);   % ��Ͳ X & Y
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_ramp_o(i,j,1),XYcoor_ramp_o(i,j,2),levelZ);   % ��Ͳ X & Y
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
            iNO,XYcoor_aisle(i,j,1),XYcoor_aisle(i,j,2),levelZ);   % ��Ͳ X & Y
    end
end
fprintf(fileID,'\n');

iNO_init2 = iNO;
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');
% �����������
% Roof_boundary = [30833,4750; 5248,41010; 45186,54754; 77970,61552; 87304,37382; 117443,7243; 114951,4750]; % ����߶�λ�㣬�����½ǵ���˳ʱ�붨λ��
levelZ = levelZaxis(end);
for i = 1:length(Roof_boundary)
    iNO = iNO+1;
    fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
        iNO,Roof_boundary(i,1),Roof_boundary(i,2),levelZ);   % ��Ͳ X & Y
end
iNO_init3 = iNO;
% ����ÿ�����
for i = 4:6 % Ŀǰȷ����¥4~6�ı���
    CoC_i = CoC(i,:);                 % centre of tower i
    levelZ_f_i = levelZ_f{i};       % facade Z axis of tower i
    facade_R_i = facade_R{i};       % facade R of tower i
    f_boundary_i = f_boundary{i};   % ��¥ i �Ӵ��ı���
    for j = 1:length(levelZ_f_i) % ����
        if facade_R_i(j) == 0   % ��������û��Ļǽ�ļ���
        else
            boundary_num = length(f_boundary_i)/4; % һ�����߻�����������
            for k = 1:boundary_num % ������
                k_s = k*4-4;
                f_b_temp1 = [f_boundary_i(k_s+1),f_boundary_i(k_s+2)];
                f_b_temp2 = [f_boundary_i(k_s+3),f_boundary_i(k_s+4)];                
                if coorPerpL(CoC_i, f_b_temp1, f_b_temp2) < facade_R_i(j) % ��������Բ�ཻ % coorPerpL�Ǵ��ߵĳ���
                    [X_temp, Y_temp, ~] = coorLxCp(CoC_i, facade_R_i(j), f_b_temp1, f_b_temp2); % �ཻ��
                    for h = 1:2
                        iNO = iNO+1;
                        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź��򣺴�0�Ƚǿ�ʼ��ʱ�룻���µ��ϡ�
                            iNO,X_temp(h),Y_temp(h),levelZ_f_i(j));   % ��Ͳ X & Y
                    end
                else
                end
            end
        end
    end
end

iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns

%% ELEMENT(frame) beams ������ҵ�����������µ�ͶӰ���Ļ���
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 1; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1���ϸֽṹQ345

% ���δ�����iPRO = 4 ������4��
fprintf(fileID,'; ���δ���\n');
ELE_iPRO = 4;
iNO = iNO_init; % ��ʼ��iNO
% �⻷��
fprintf(fileID,';   �µ�����\n');
for k = 1:2 % ��ҵ�㡢�����
    for i = 1:lengthXYcoor_ramp % �µ�ͶӰһȦ�Ľڵ���
        for j = 1:circle_num % �µ���/�⻷��/�ߵ��⻷��
            iEL = iEL+1;
            iN1 = iNO+j+(i-1)*circle_num+(k-1)*lengthXYcoor_ramp*circle_num;	% ��Ͳ�ϵĵ�/�µ����ϵĵ�
            if i == lengthXYcoor_ramp
                iN2 = iN1+circle_num-lengthXYcoor_ramp*circle_num;
            else
                iN2 = iN1+circle_num;    % �µ����ϵĵ�/�µ���/�ߵ����ϵĵ�
            end
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % ����Ԫ�������ڵ��
                ELE_ANGLE, ELE_iSUB);
        end
    end
end
fprintf(fileID,'\n');

fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');
% ���滷�α�����iPRO = 4 ������4��
fprintf(fileID,'; ���滷�α���\n');
ELE_iPRO = 4;
iNO = iNO_init2; % ��ʼ��iNO
for i = 1:length(Roof_boundary)
    iEL = iEL+1;
    iN1 = iNO+i;	% һ�ֳ�����
    if i == length(Roof_boundary)
        iN2 = iNO+1;
    else
        iN2 = iN1+1;    %
    end
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % ����Ԫ�������ڵ��
        ELE_ANGLE, ELE_iSUB);
end
fprintf(fileID,'\n');
% ÿ�㻷�α�����iPRO = 4 ������4��
fprintf(fileID,'; ÿ�㻷�α���\n');
ELE_iPRO = 4;
iNO = iNO_init3; % ��ʼ��iNO
iN1 = iNO-1;
while iN1 < iNO_end-1
    iEL = iEL+1;
    iN1 = iN1+2;	% 12,34,56...
    iN2 = iN1+1;
    fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
        iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
        iN1, iN2,...    % ����Ԫ�������ڵ��
        ELE_ANGLE, ELE_iSUB);
end
fprintf(fileID,'\n');

iEL_end = iEL;
end
