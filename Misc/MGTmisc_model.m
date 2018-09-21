%% function
% MGTmisc model
% 
% Xu Yi, 2018

%%
function MGTmisc_model(fileID)
%% append TOWERs
% initial conditions
iNO = 15000;    % �ڵ�ų�ʼ�� ԭģ�ͼ�����14398 ��15000��ʼ
iEL = 0;
% car_num = 8;    % Բ��ÿ��ͣ����
tower_colu_num = 4; % ͣ��С����������
% 10����¥��Բ�� % centre of a circle
CoC_towerC1   = [57400,36000];	% centre parking tower
CoC_towerS2   = [26225,27850];	% side parking tower 1
CoC_towerS3   = [86800,12000];	% side parking tower 2
CoC_elevator4 = [61600,10800];	% 1 elevator
CoC_stair5    = [36400,42300];	% stairs
CoC_stair6    = [78400,31800];	% stairs
CoC_side7     = [75366,57542];	% side type 1
CoC_side8     = [117104,5443];	% side type 2 / 2 columns in the edge
CoC_side9     = [30655, 5307];	% side type 3 / a column in the centre
CoC_side10    = [ 5514,38672];	% side type 1

Edge_side9 = [29407,2900]; Edge_side10 = [3769,39236]; Edge_side7 = [76602,59809];
Edge_side8 = [116115,5443; 113571,2900]; Edge_North = [43773,53002]; Edge_East = [85996,35635];
Roof_boundary = [Edge_side9; Edge_side10; Edge_North; Edge_side7; Edge_East; Edge_side8]; % ����߶�λ�㣬�����½ǵ���˳ʱ�붨λ��

Deg_ref = atan( (CoC_elevator4(1)-CoC_towerC1(1))/(CoC_elevator4(2)-CoC_towerC1(2)) );	% �ɵ���Բ��������Բ������ȷ��
% Deg_towerC1 = 0; % degree of centre tower ���ǵ���һ�㳵�������⣬������ת��
Deg_towerS2 = 0; % degree of side tower 1 �����
Deg_towerS3 = 0; % degree of side tower 2 �ұ���
Deg_elevator4 = Deg_ref;	% degree of elevator4 �ɵ���Բ��������Բ������ȷ��
Deg_stair5 = pi/4 + Deg_ref;   % �ɵ���Բ��������Բ������ȷ��
Deg_stair6 = -acot( (CoC_stair6(1)-CoC_towerC1(1))/(CoC_stair6(2)-CoC_towerC1(2)) )+pi;	% degree of stair 6 ��¥��6Բ��������Բ������ȷ��
% Deg_side7 = pi/4 + Deg_towerC1;    % �������Ƕ�ȷ�� % Deg_side8 = -pi/6;    % ���� % Deg_side9 = -pi/2 + Deg_towerS2;   % ��S2���Ƕ�ȷ�� % Deg_side10 = -pi/6;    % ����

% facade_tower1_R = [zeros(3,1); 10298; 9802; 9700; 10032; 10830; 12105; 14256; 18500]; % Ļǽ������500 Ϊ�ṹ����
facade_tower2_R = [9626, 11247];
facade_tower3_R = [9274, 10734];
facade_ele4_R = [7816, 9202];
facade_stair5_R = [8914, 10255];
facade_stair6_R = [7954, 9264];
facade_side7_R = [7308, 9510];
facade_side8_R = [9597, 11684];
facade_side9_R = [7328, 9683];
facade_side10_R = [6566, 8946];

towerS_column_coor = [ 3125,4150; 3850,3500 ]; % ����ͣ��С��4�������γɵľ��εĳ��Ϳ�
levelaxis = [18700, 20200];
levelPstart = [1, 1, 1, 1, 1]; % ͣ����¥�㣬��¥�����йء���һ��Ϊ��Ͳ���ڶ���ΪСͲ(����û��)�� �������ĸ�Ϊ¥��ר��(3/��ƽ̨��4/����ƽ̨)�������Ϊ����ר�á�

stairColu_num = 4;  % ¥����Ͳ������
stairW = 3950; % ¥�ݿ�
stairL = 4300; % ¥�ݳ�������̤��ǰ������
elevatorColu_num = 8;  % ����Ͳ����Ͳ������ (�����м仹��һ���ڵ㣬��7+1=8)
elevatorR = 4250; % ����Ͳ�İ뾶

Arc_itvl = 1000; % ���塰��ֱ�����������ֱ�߶γ��ȡ�
CAR = 0;
OFFICE = 0;
ROOF = 0;

%% 10 towers
% fprintf(fileID,'; ��1\n');
% iNO_towerC1_init = iNO;
% [iNO, iEL] = MGT_facade_C1(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, facade_tower1_R, levelTaxis, levelPstart(1), iNO_towerC1_init, Arc_itvl);

fprintf(fileID,'; ��2\n');
iNO_towerS2_init = iNO;
[iNO, iEL] = MGTmisc_facade_S2(fileID, iNO, iEL, tower_colu_num, CoC_towerS2, Deg_towerS2, towerS_column_coor(1,:), facade_tower2_R, levelaxis, levelPstart, iNO_towerS2_init, Arc_itvl);
fprintf(fileID,'; ��3\n');
iNO_towerS3_init = iNO;
[iNO, iEL] = MGTmisc_facade_S3(fileID, iNO, iEL, tower_colu_num, CoC_towerS3, Deg_towerS3, towerS_column_coor(2,:), facade_tower3_R, levelaxis, levelPstart, iNO_towerS3_init, Arc_itvl);

fprintf(fileID,'; ��4\n');
[iNO, iEL] = MGTmisc_elevator(fileID, iNO, iEL, CoC_elevator4, Deg_elevator4, facade_ele4_R, levelaxis, levelPstart(4), elevatorColu_num, elevatorR, CAR, OFFICE, ROOF, Arc_itvl);

fprintf(fileID,'; ��5\n');
[iNO, iEL] = MGTmisc_stair(fileID, iNO, iEL, CoC_stair5, Deg_stair5, facade_stair5_R, levelaxis, levelPstart(4), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF, Arc_itvl);
fprintf(fileID,'; ��6\n');
[iNO, iEL] = MGTmisc_stair(fileID, iNO, iEL, CoC_stair6, Deg_stair6, facade_stair6_R, levelaxis, levelPstart(4), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF, Arc_itvl);

fprintf(fileID,'; ��7\n');
[iNO, iEL] = MGTmisc_side(fileID, iNO, iEL, CoC_side7, facade_side7_R, levelaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF, 7, Arc_itvl);

fprintf(fileID,'; ��8\n');
[iNO, iEL] = MGTmisc_side8(fileID, iNO, iEL, CoC_side8, facade_side8_R, levelaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF, Arc_itvl);
fprintf(fileID,'; ��9\n');
[iNO, iEL] = MGTmisc_side(fileID, iNO, iEL, CoC_side9, facade_side9_R, levelaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF, 9, Arc_itvl);

fprintf(fileID,'; ��10\n');
[~, ~] = MGTmisc_side(fileID, iNO, iEL, CoC_side10, facade_side10_R, levelaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF, 10, Arc_itvl);

%%
fprintf(fileID,'\n');

end