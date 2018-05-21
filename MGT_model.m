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
iNO = 0;    % 节点号初始化
iEL = 0;    % 单元号初始化
car_num = 8;    % 圆塔每层停车数
% 10个塔楼的圆心 % centre of a circle
CoC_towerC1 = [58800,37800];	% centre parking tower
CoC_towerS2 = [27300,29400];	% side parking tower 1
CoC_towerS3 = [88200,13800];	% side parking tower 2

CoC_elevator4 = [63000,12600];	% 1 elevator

CoC_stair5 = [36977,44429];	% stairs
CoC_stair6 = [79800,33600];	% stairs

CoC_side7 = [76552,58934];	% side type 1
CoC_side10 = [6914,40472];	% side type 1

CoC_side8 = [118504,7243];	% side type 2 / 2 columns in the edge
CoC_side9 = [32055,7107];	% side type 3 / a column in the centre

Edge_side9 = [30833,4750]; Edge_side10 = [5248,41010]; Edge_side7 = [77970,61552];
Edge_side8 = [117443,7243; 114951,4750]; Edge_North = [45186,54754]; Edge_East = [87304,37382];
Roof_boundary = [Edge_side9; Edge_side10; Edge_North; Edge_side7; Edge_East; Edge_side8]; % 外边线定位点，从左下角点起，顺时针定位点

Deg_towerC1 = atan( (CoC_elevator4(1)-CoC_towerC1(1))/(CoC_elevator4(2)-CoC_towerC1(2)) );	% degree of centre tower 由电梯圆心与主塔圆心连线确定
Deg_towerS2 = -acot( (CoC_stair5(1)-CoC_towerS2(1))/(CoC_stair5(2)-CoC_towerS2(2)) ) + pi/4;	% degree of side tower 1 由中上楼梯圆心与左边塔圆心连线确定
Deg_towerS3 = 0; % degree of side tower 2 右边塔
Deg_elevator4 = Deg_towerC1;    % 由主塔角度确定
Deg_stair5 = pi/4 + Deg_towerS2;   % 由S2塔角度确定
Deg_stair6 = -acot( (CoC_stair6(1)-CoC_towerC1(1))/(CoC_stair6(2)-CoC_towerC1(2)) );	% degree of stair 6 由楼梯6圆心与主塔圆心连线确定
% Deg_side7 = pi/4 + Deg_towerC1;    % 由主塔角度确定
Deg_side8 = -pi/6;    % 待定
Deg_side9 = -pi/2 + Deg_towerS2;   % 由S2塔角度确定
Deg_side10 = -pi/6;    % 待定

facade_side7_R = [zeros(6,1); 6820; 5839; 5162; 4734; 4529; 4534; 4750; 5189; 5880; 6878; 8473; 11040; 16200];
facade_side10_R = [zeros(6,1); 6927; 5477; 4483; 3827; 3452; 3333; 3460; 3843; 4509; 5514; 7163; 9860; 15333];

tube_innerR = 3950;
tube_outerR = 8500;
levelTaxis = [-6200:2450:3600, 5800:2200:12400, 15000:2600:17600, 19800, 22350];  % 塔楼楼层标高 原-2100:2100:23100
levelSaxis = [-6200, -4524, -3016, -1508, 0:1600:17600, 19360, 21120, 22350];  % 楼梯楼层标高
levelPstart = [5, length(levelTaxis(:))-2, 7]; % 停车的楼层，与楼层数有关。第一个为大筒，第二个为小筒。 第三个为楼梯专用。

stairColu_num = 4;  % 楼梯内筒柱数量
stairL = 3300; % 楼梯长，即沿踏步前进方向长
stairW = 3950; % 楼梯宽
stairB = 1500; % 楼梯梯板宽(暂定)
elevatorColu_num = 8;  % 电梯筒的内筒柱数量 (电梯中间还有一个节点，故7+1=8)

%% 10 towers
iNO_towerC1_init = iNO;
[iNO, iEL] = MGT_tower(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, tube_outerR, levelTaxis, levelPstart(1), CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_ramp(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, tube_outerR, levelTaxis, levelPstart, iNO_towerC1_init);
[iNO, iEL] = MGT_Misc(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, tube_outerR, levelTaxis, levelPstart, Roof_boundary, iNO_towerC1_init);

iNO_towerS2_init = iNO;
[iNO, iEL] = MGT_tower(fileID, iNO, iEL, car_num, CoC_towerS2, Deg_towerS2, tube_innerR, tube_outerR, levelTaxis, levelPstart(2), CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_facade_tower(fileID, iNO, iEL, car_num, CoC_towerS2, Deg_towerS2, tube_innerR, levelTaxis, levelPstart, iNO_towerS2_init);

iNO_towerS3_init = iNO;
[iNO, iEL] = MGT_tower(fileID, iNO, iEL, car_num, CoC_towerS3, Deg_towerS3, tube_innerR, tube_outerR, levelTaxis, levelPstart(2), CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_facade_tower(fileID, iNO, iEL, car_num, CoC_towerS3, Deg_towerS3, tube_innerR, levelTaxis, levelPstart, iNO_towerS3_init);

[iNO, iEL] = MGT_elevator(fileID, iNO, iEL, CoC_elevator4, Deg_elevator4, levelSaxis, levelPstart(3), elevatorColu_num, CAR, OFFICE, ROOF);

iNO_stair5_init = iNO;
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_stair5, Deg_stair5, levelSaxis, levelPstart(3), stairColu_num, stairL, stairW, stairB, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_facade_stair(fileID, iNO, iEL, stairColu_num, CoC_stair5, Deg_stair5, stairL, stairW, levelSaxis, levelPstart(3), iNO_stair5_init, 5);
iNO_stair6_init = iNO;
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_stair6, Deg_stair6, levelSaxis, levelPstart(3), stairColu_num, stairL, stairW, stairB, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_facade_stair(fileID, iNO, iEL, stairColu_num, CoC_stair6, Deg_stair6, stairL, stairW, levelSaxis, levelPstart(3), iNO_stair6_init, 6);

[iNO, iEL] = MGT_side(fileID, iNO, iEL, CoC_side7, facade_side7_R, levelSaxis, levelPstart(3), Roof_boundary, CAR, OFFICE, ROOF, 7);
[iNO, iEL] = MGT_side(fileID, iNO, iEL, CoC_side10, facade_side10_R, levelSaxis, levelPstart(3), Roof_boundary, CAR, OFFICE, ROOF, 10);

iNO_stair8_init = iNO;
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_side8, Deg_side8, levelSaxis, levelPstart(3), stairColu_num, stairL, stairW, stairB, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_facade_stair(fileID, iNO, iEL, stairColu_num, CoC_side8, Deg_side8, stairL, stairW, levelSaxis, levelPstart(3), iNO_stair8_init, 8);
iNO_stair9_init = iNO;
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_side9, Deg_side9, levelSaxis, levelPstart(3), stairColu_num, stairL, stairW, stairB, CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_facade_stair(fileID, iNO, iEL, stairColu_num, CoC_side9, Deg_side9, stairL, stairW, levelSaxis, levelPstart(3), iNO_stair9_init, 9);

%%


%%
fprintf(fileID,'\n');

end