%% function
% MGT boundary beams.
%
% Xu Yi, 2018

%%
function [iNO_end, iEL_end] = MGT_boundary(fileID, iNO, iEL, ~, CoC_tower, ~, ~, ~, levelZaxis, ~, Roof_boundary, ~,...
    CoC_towerS2, CoC_towerS3, CoC_elevator4, CoC_stair5, CoC_stair6, facade_tower2_R, facade_tower3_R, facade_ele4_R, facade_stair5_R, facade_stair6_R,...
    levelTaxis, levelTSaxis, levelSaxis_f)    % ע�������levelPstart��1x2����
%%
CoC = [CoC_tower; CoC_towerS2; CoC_towerS3; CoC_elevator4; CoC_stair5; CoC_stair6];                     % centre of tower1~6
facade_R = {[0,0]; facade_tower2_R; facade_tower3_R; facade_ele4_R; facade_stair5_R; facade_stair6_R};    % facade R of tower1~6
levelZ_f = {levelTaxis; levelTSaxis; levelTSaxis; levelSaxis_f; levelSaxis_f; levelSaxis_f};            % facade Z axis of tower1~6
f_boundary = {[0,0,0,0]; [Roof_boundary(1,:),Roof_boundary(2,:)]; [Roof_boundary(7,:),Roof_boundary(1,:)]; [Roof_boundary(7,:),Roof_boundary(1,:)];...
    [Roof_boundary(2,:),Roof_boundary(3,:),Roof_boundary(3,:),Roof_boundary(4,:)]; [Roof_boundary(4,:),Roof_boundary(5,:),Roof_boundary(5,:),Roof_boundary(6,:)]};  % 1~6��¥�Ӵ��ı���
%% NODE ������ҵ�����������µ�ͶӰ���Ľڵ�
iNO_init1 = iNO;
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
iNO_init2 = iNO;
% ����ÿ�����
for i = 2:6 % Ŀǰȷ����¥2~6�ı���
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

% ���滷�α�����iPRO = 4 ������4��
fprintf(fileID,'; ���滷�α���\n');
ELE_iPRO = 4;
iNO = iNO_init1; % ��ʼ��iNO
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
iNO = iNO_init2; % ��ʼ��iNO
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
