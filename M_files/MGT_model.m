%% function
% MGT model
% 
% Xu Yi, 23rd March 2018
% Xu Yi, 22nd April 2018, revised

%%
function MGT_model(fileID)
%% SECTION
SEC_P450x12_num = '1'; SEC_P500x14_num = '2'; SEC_B500X300X12_num = '3'; SEC_B350x250x8_num = '4'; SEC_P245x12_num = '5';
SEC_P450x12_SNAME = 'P 450x12'; SEC_P500x14_SNAME = 'P 500x14'; SEC_B500X300X12_SNAME = 'B 500x300x12'; SEC_B350x250x8_SNAME = 'B 350x250x8'; SEC_P245x12_SNAME = 'P 245x12';
SEC_TYPE = 'DBUSER'; SEC_OFFSET = 'CC'; SEC_0 = '0'; SEC_1 = '1'; SEC_YES = 'YES'; SEC_NO = 'NO'; SEC_P = 'P';  SEC_B = 'B';  SEC_GB = 'GB-YB05';

fprintf(fileID,'*SECTION    ; Section\n');
fprintf(fileID,'; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, [DATA1], [DATA2]                    ; 1st line - DB/USER\n; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, BLT, D1, ..., D8, iCEL              ; 1st line - VALUE\n;       AREA, ASy, ASz, Ixx, Iyy, Izz                                               ; 2nd line\n;       CyP, CyM, CzP, CzM, QyB, QzB, PERI_OUT, PERI_IN, Cy, Cz                     ; 3rd line\n;       Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, Zyy, Zzz                                    ; 4th line\n; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, ELAST, DEN, POIS, POIC, SF, THERMAL ; 1st line - SRC\n;       D1, D2, [SRC]                                                               ; 2nd line\n; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, 1, DB, NAME1, NAME2, D1, D2         ; 1st line - COMBINED\n; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, 2, D11, D12, D13, D14, D15, D21, D22, D23, D24\n; iSEC, TYPE, SNAME, [OFFSET2], bSD, bWE, SHAPE, iyVAR, izVAR, STYPE                ; 1st line - TAPERED\n;       DB, NAME1, NAME2                                                            ; 2nd line(STYPE=DB)\n;       [DIM1], [DIM2]                                                              ; 2nd line(STYPE=USER)\n;       D11, D12, D13, D14, D15, D16, D17, D18                                      ; 2nd line(STYPE=VALUE)\n;       AREA1, ASy1, ASz1, Ixx1, Iyy1, Izz1                                         ; 3rd line(STYPE=VALUE)\n;       CyP1, CyM1, CzP1, CzM1, QyB1, QzB1, PERI_OUT1, PERI_IN1, Cy1, Cz1           ; 4th line(STYPE=VALUE)\n;       Y11, Y12, Y13, Y14, Z11, Z12, Z13, Z14, Zyy1, Zyy2                          ; 5th line(STYPE=VALUE)\n;       D21, D22, D23, D24, D25, D26, D27, D28                                      ; 6th line(STYPE=VALUE)\n;       AREA2, ASy2, ASz2, Ixx2, Iyy2, Izz2                                         ; 7th line(STYPE=VALUE)\n;       CyP2, CyM2, CzP2, CzM2, QyB2, QzB2, PERI_OUT2, PERI_IN2, Cy2, Cz2           ; 8th line(STYPE=VALUE)\n;       Y21, Y22, Y23, Y24, Z21, Z22, Z23, Z24, Zyy2, Zzz2                          ; 9th line(STYPE=VALUE)\n; [DATA1] : 1, DB, NAME or 2, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10\n; [DATA2] : CCSHAPE or iCEL or iN1, iN2\n; [SRC]  : 1, DB, NAME1, NAME2 or 2, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, iN1, iN2\n; [DIM1], [DIM2] : D1, D2, D3, D4, D5, D6, D7, D8\n; [OFFSET] : OFFSET, iCENT, iREF, iHORZ, HUSER, iVERT, VUSER\n; [OFFSET2]: OFFSET, iCENT, iREF, iHORZ, HUSERI, HUSERJ, iVERT, VUSERI, VUSERJ\n');

fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
    SEC_P450x12_num,SEC_TYPE,SEC_P450x12_SNAME,...
    SEC_OFFSET,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_YES,SEC_NO,...
    SEC_P,SEC_1,SEC_GB,SEC_P450x12_SNAME);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
    SEC_P500x14_num,SEC_TYPE,SEC_P500x14_SNAME,...
    SEC_OFFSET,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_YES,SEC_NO,...
    SEC_P,SEC_1,SEC_GB,SEC_P500x14_SNAME);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
    SEC_B500X300X12_num,SEC_TYPE,SEC_B500X300X12_SNAME,...
    SEC_OFFSET,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_YES,SEC_NO,...
    SEC_B,SEC_1,SEC_GB,SEC_B500X300X12_SNAME);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
    SEC_B350x250x8_num,SEC_TYPE,SEC_B350x250x8_SNAME,...
    SEC_OFFSET,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_YES,SEC_NO,...
    SEC_B,SEC_1,SEC_GB,SEC_B350x250x8_SNAME);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
    SEC_P245x12_num,SEC_TYPE,SEC_P245x12_SNAME,...
    SEC_OFFSET,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_YES,SEC_NO,...
    SEC_P,SEC_1,SEC_GB,SEC_P245x12_SNAME);
fprintf(fileID,'\n');

%% THICKNESS
THICK_120_num = '1'; THICK_1_num = '2';
THICK_120 = '120'; THICK_1 = '1';
THICK_TYPE = 'VALUE'; THICK_bSAME = 'YES';
THICK_bOFFSET = 'NO'; THICK_OFFTYPE = '0'; THICK_VALUE = '0';

fprintf(fileID,'*THICKNESS    ; Thickness\n');
fprintf(fileID,'; iTHK, TYPE, bSAME, THIK-IN, THIK-OUT, bOFFSET, OFFTYPE, VALUE ; TYPE=VALUE\n; iTHK, TYPE, SUBTYPE, RPOS, WEIGHT                             ; TYPE=STIFFENED, SUBTYPE=VALUE\n;       SHAPE, THIK-IN, THIK-OUT, HU, HL                        ;      for yz section\n;       SHAPE, THIK-IN, THIK-OUT, HU, HL                        ;      for xz section\n; iTHK, TYPE, SUBTYPE, RPOS, PLATETHIK                          ; TYPE=STIFFENED, SUBTYPE=USER\n;       bRIB {, SHAPE, DIST, SIZE1, SIZE2, ..., SIZE6}          ;      for yz section\n;       bRIB {, SHAPE, DIST, SIZE2, SIZE2, ..., SIZE6}          ;      for xz section\n; iTHK, TYPE, SUBTYPE, RPOS, PLATETHIK, DBNAME                  ; TYPE=STIFFENED, SUBTYPE=DB\n;       bRIB {, SHAPE, DIST, SNAME}                             ;      for yz section\n;       bRIB {, SHAPE, DIST, SNAME}                             ;      for xz section\n');
fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s\n',...
    THICK_120_num,THICK_TYPE,THICK_bSAME,THICK_120,THICK_120,THICK_bOFFSET,THICK_OFFTYPE,THICK_VALUE);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s\n',...
    THICK_1_num,THICK_TYPE,THICK_bSAME,THICK_1,THICK_1,THICK_bOFFSET,THICK_OFFTYPE,THICK_VALUE);
fprintf(fileID,'\n');

%% FLOADTYPE
LCNAME_DL = 'DL'; LCNAME_LL = 'LL';
bSBU_YES = 'YES'; bSBU_NO = 'NO';

CAR = 'CAR'; OFFICE = 'OFFICE'; ROOF = 'ROOF';
CAR_DL = '-2e-006'; CAR_LL = '-4e-006';
OFFICE_DL = '-1.5e-006'; OFFICE_LL = '-3.5e-006';
ROOF_DL = '-4e-006'; ROOF_LL = '-3.5e-006';

fprintf(fileID,'*FLOADTYPE    ; Define Floor Load Type\n');
fprintf(fileID,'; NAME, DESC                                           ; 1st line\n; LCNAME1, FLOAD1, bSBU1, ..., LCNAME8, FLOAD8, bSBU8  ; 2nd line\n');
fprintf(fileID,'   %s, \n', CAR);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s\n',...
                LCNAME_DL, CAR_DL, bSBU_YES, LCNAME_LL, CAR_LL, bSBU_NO);
fprintf(fileID,'   %s, \n', OFFICE);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s\n',...
                LCNAME_DL, OFFICE_DL, bSBU_YES, LCNAME_LL, OFFICE_LL, bSBU_NO);
fprintf(fileID,'   %s, \n', ROOF);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s\n',...
                LCNAME_DL, ROOF_DL, bSBU_YES, LCNAME_LL, ROOF_LL, bSBU_NO);
fprintf(fileID,'\n');

%% append TOWERs
% initial conditions
iNO = 0;    % �ڵ�ų�ʼ��
iEL = 0;    % ��Ԫ�ų�ʼ��
car_num = 8;    % Բ��ÿ��ͣ����
% 10����¥��Բ�� % centre of a circle
CoC_towerC1 = [58800,37800]; % centre parking tower
CoC_towerS2 = [27300,29400]; % side parking tower 1
CoC_towerS3 = [88200,13800]; % side parking tower 2

CoC_stairs5 = [37000,44400]; % stairs
CoC_stairs6 = [79800,33600]; % stairs
CoC_stairs7 = [76550,58935]; % stairs ����
CoC_stairs9 = [32050,7450]; % stairs ����
            
CoC_elevator4 = [63000,12600]; % 1 elevator
CoC_side8 = [113970,8165]; % side
CoC_side10 = [8400,40000]; % side

Deg_towerC1 = atan( (CoC_elevator4(1)-CoC_towerC1(1))/(CoC_elevator4(2)-CoC_towerC1(2)) );	% degree of centre tower �ɵ���Բ��������Բ������ȷ��
Deg_towerS2 = -acot( (CoC_stairs5(1)-CoC_towerS2(1))/(CoC_stairs5(2)-CoC_towerS2(2)) ) + pi/4;	% degree of side tower 1 ������¥��Բ���������Բ������ȷ��
Deg_towerS3 = 0; % degree of side tower 2 �ұ���
Deg_stair5 = pi/4 + Deg_towerS2;   % ��S2���Ƕ�ȷ��
Deg_stair6 = -acot( (CoC_stairs6(1)-CoC_towerC1(1))/(CoC_stairs6(2)-CoC_towerC1(2)) );	% degree of stair 6 ��¥��6Բ��������Բ������ȷ��
Deg_stair7 = pi/4 + Deg_towerC1;    % �������Ƕ�ȷ��
Deg_stair9 = -pi/2 + Deg_towerS2;   % ��S2���Ƕ�ȷ��
Deg_elevator4 = Deg_towerC1;    % �������Ƕ�ȷ��
Deg_stair10 = -pi/6;    % ����
Deg_stair8 = -pi/6;    % ����

tube_innerR = 3950;
tube_outerR = 8500;
levelZaxis = -2100:2100:23100;  % ¥����
levelPstart = [3,length(levelZaxis(:))-2]; % ͣ����¥�㣬��¥�����йء���һ��Ϊ��Ͳ��������ΪСͲ��

stairColu_num = 4;  % ¥����Ͳ������
stairL = 2750; % ¥�ݳ�������̤��ǰ������
stairW = 3150; % ¥�ݿ�
elevatorColu_num = 7;  % ����Ͳ����Ͳ������

%% 10 towers
iNO_towerC_init = iNO;
[iNO, iEL] = MGT_tower(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, tube_outerR, levelZaxis, levelPstart(1), CAR, OFFICE, ROOF);
iNO_towerS1_init = iNO;
[iNO, iEL] = MGT_tower(fileID, iNO, iEL, car_num, CoC_towerS2, Deg_towerS2, tube_innerR, tube_outerR, levelZaxis, levelPstart(2), CAR, OFFICE, ROOF);
iNO_towerS2_init = iNO;
[iNO, iEL] = MGT_tower(fileID, iNO, iEL, car_num, CoC_towerS3, Deg_towerS3, tube_innerR, tube_outerR, levelZaxis, levelPstart(2), CAR, OFFICE, ROOF);

[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_stairs5, Deg_stair5, levelZaxis, levelPstart(1), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_stairs6, Deg_stair6, levelZaxis, levelPstart(1), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_stairs7, Deg_stair7, levelZaxis, levelPstart(1), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_stairs9, Deg_stair9, levelZaxis, levelPstart(1), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF);

[iNO, iEL] = MGT_elevator(fileID, iNO, iEL, CoC_elevator4, Deg_elevator4, levelZaxis, levelPstart(1), elevatorColu_num, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_side10, Deg_stair10, levelZaxis, levelPstart(1), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_side8, Deg_stair8, levelZaxis, levelPstart(1), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF);

%%


%%
fprintf(fileID,'\n');

end