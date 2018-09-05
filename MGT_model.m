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
tower_colu_num = 4; % 停车小塔的柱子数
% 10个塔楼的圆心 % centre of a circle
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
Roof_boundary = [Edge_side9; Edge_side10; Edge_North; Edge_side7; Edge_East; Edge_side8]; % 外边线定位点，从左下角点起，顺时针定位点

Deg_ref = atan( (CoC_elevator4(1)-CoC_towerC1(1))/(CoC_elevator4(2)-CoC_towerC1(2)) );	% 由电梯圆心与主塔圆心连线确定
Deg_towerC1 = 0; % degree of centre tower 考虑地下一层车辆进塔库，该塔不转。
Deg_towerS2 = 0; % degree of side tower 1 左边塔
Deg_towerS3 = 0; % degree of side tower 2 右边塔
Deg_elevator4 = Deg_ref;	% degree of elevator4 由电梯圆心与主塔圆心连线确定
Deg_stair5 = pi/4 + Deg_ref;   % 由电梯圆心与主塔圆心连线确定
Deg_stair6 = -acot( (CoC_stair6(1)-CoC_towerC1(1))/(CoC_stair6(2)-CoC_towerC1(2)) )+pi;	% degree of stair 6 由楼梯6圆心与主塔圆心连线确定
% Deg_side7 = pi/4 + Deg_towerC1;    % 由主塔角度确定 % Deg_side8 = -pi/6;    % 待定 % Deg_side9 = -pi/2 + Deg_towerS2;   % 由S2塔角度确定 % Deg_side10 = -pi/6;    % 待定

facade_tower1_R = [zeros(3,1); 10298; 9802; 9700; 10032; 10830; 12105; 14346; 18500]; % 幕墙线内退500 为结构边线
% 备份 facade_tower2_R = [zeros(5,1); 5328; 4597; 4671; 5582; 7205; 8970; 12715];
facade_tower2_R = [zeros(3,1); 7220; 6352; 5741; 5807; 6756; 8470; 12715];
facade_tower3_R = [zeros(3,1); 7220; 6352; 5741; 5807; 6756; 8145; 12030];
% 备份 facade_ele4_R = [zeros(6,1); 8273; 7304; 6661; 6286; 6155; 6260; 6610; 7233; 8188; 9610; 12037; 15239; 17512];
% 备份 facade_ele4_R = [zeros(3,1); 5411; 3982; 3539; 4064; 5580; 6705; 10044; 10400];
facade_ele4_R = [zeros(3,1); 5411; 4500; 4500; 4500; 5580; 6705; 10400];
% 备份 facade_stair5_R = [zeros(6,1); 8281; 7316; 6675; 6303; 6172; 6278; 6627; 7249; 8202; 9619; 12004; 15218; 17618];
facade_stair5_R = [zeros(3,1); 6960; 5594; 5088; 5530; 6918; 7962; 11415];
facade_stair6_R = [zeros(3,1); 5854; 4683; 4244; 4705; 6054; 7059; 10370];
facade_side7_R = [zeros(3,1); 8544; 5134; 3939; 3636; 4357; 5875; 12320];
facade_side8_R = [zeros(3,1); 9597; 7137; 6018; 5889; 6714; 8227; 14060];
facade_side9_R = [zeros(3,1); 7264; 4797; 3688; 3532; 4323; 5862; 12500];
facade_side10_R = [zeros(3,1); 6433; 3986; 2966; 2870; 3645; 5130; 11845];

tube_innerR = 3950;
tube_outerR = 8500;
towerS_column_coor = [ 3125,4150; 3850,3500 ]; % 两个停车小塔4个柱子形成的矩形的长和宽
levelTaxis = [-6300:3100:-100, 2300:2400:14300, 17300, 22000];  % 塔楼楼层标高 -6300(-1)/3500/17500(commercial)/22000(Roof)
levelTSaxis = [-6300:3100:-100, 2900:3000:14900, 17300, 22000]; % 
% levelTSaxis = [-6300, 0:3500:17500, 22000];  % 小塔楼楼层标高 -6300(-1)/3500/17500(commercial)/22000(Roof)
levelSaxis = [-6300, -4624, -3116, -1608, -100:1600:15900, 17300, 19260, 21020, 22000];  % 楼梯楼层标高
levelSaxis_f = [-6300, -3116, -100:3200:15900, 17300, 22000];  % 楼梯幕墙支撑梁标高
levelSDaxis = levelTSaxis;  % 角塔楼层标高
levelPstart = [4, length(levelTaxis(:))-2, 7, 4, 4]; % 停车的楼层，与楼层数有关。第一个为大筒，第二个为小筒(现在没用)。 第三、四个为楼梯专用(3/含平台，4/不含平台)。第五个为角塔专用。

stairColu_num = 4;  % 楼梯内筒柱数量
stairW = 3950; % 楼梯宽
stairL = 4300; % 楼梯长，即沿踏步前进方向长
stairB = 1551; % 楼梯梯板宽(暂定) 1975*pi/4 半圆与矩形的等面积代换
elevatorColu_num = 8;  % 电梯筒的内筒柱数量 (电梯中间还有一个节点，故7+1=8)
elevatorR = 4250; % 电梯筒的半径

Arc_itvl = 1000; % 定义“以直代曲”的最大直线段长度。

%% 10 towers
fprintf(fileID,'; 塔1\n');
iNO_towerC1_init = iNO;
[iNO, iEL] = MGT_tower(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, tube_outerR, levelTaxis, levelPstart(1), CAR, OFFICE, ROOF);
% [iNO, iEL] = MGT_ramp(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, tube_outerR, levelTaxis, levelPstart, iNO_towerC1_init);
[iNO, iEL] = MGT_facade_C1(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, facade_tower1_R, levelTaxis, levelPstart(1), iNO_towerC1_init, Arc_itvl);
[iNO, iEL] = MGT_boundary(fileID, iNO, iEL, car_num, CoC_towerC1, Deg_towerC1, tube_innerR, tube_outerR, levelTaxis, levelPstart, Roof_boundary, iNO_towerC1_init, ...
                CoC_towerS2, CoC_towerS3, CoC_elevator4, CoC_stair5, CoC_stair6,...
                facade_tower2_R, facade_tower3_R, facade_ele4_R, facade_stair5_R, facade_stair6_R,...
                levelTaxis, levelTSaxis, levelSaxis_f);

fprintf(fileID,'; 塔2\n');
iNO_towerS2_init = iNO;
[iNO, iEL] = MGT_tower_S(fileID, iNO, iEL, tower_colu_num, CoC_towerS2, Deg_towerS2, towerS_column_coor(1,:), levelTSaxis, levelPstart(2), CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_facade_S2(fileID, iNO, iEL, tower_colu_num, CoC_towerS2, Deg_towerS2, towerS_column_coor(1,:), facade_tower2_R, levelTSaxis, levelPstart, iNO_towerS2_init, Arc_itvl);
fprintf(fileID,'; 塔3\n');
iNO_towerS3_init = iNO;
[iNO, iEL] = MGT_tower_S(fileID, iNO, iEL, tower_colu_num, CoC_towerS3, Deg_towerS3, towerS_column_coor(2,:), levelTSaxis, levelPstart(2), CAR, OFFICE, ROOF);
[iNO, iEL] = MGT_facade_S3(fileID, iNO, iEL, tower_colu_num, CoC_towerS3, Deg_towerS3, towerS_column_coor(2,:), facade_tower3_R, levelTSaxis, levelPstart, iNO_towerS3_init, Arc_itvl);

fprintf(fileID,'; 塔4\n');
[iNO, iEL] = MGT_elevator(fileID, iNO, iEL, CoC_elevator4, Deg_elevator4, facade_ele4_R, levelSaxis_f, levelPstart(4), elevatorColu_num, elevatorR, CAR, OFFICE, ROOF, Arc_itvl);
[iNO, iEL] = MGT_elevator_platform(fileID, iNO, iEL, CoC_elevator4, Deg_elevator4, levelSaxis, elevatorColu_num, elevatorR, CAR, OFFICE, ROOF);

fprintf(fileID,'; 塔5\n');
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_stair5, Deg_stair5, facade_stair5_R, levelSaxis_f, levelPstart(4), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF, Arc_itvl);
[iNO, iEL] = MGT_stair_platform(fileID, iNO, iEL, CoC_stair5, Deg_stair5, levelSaxis, stairColu_num, stairL, stairW, stairB, CAR, OFFICE, ROOF);
fprintf(fileID,'; 塔6\n');
[iNO, iEL] = MGT_stair(fileID, iNO, iEL, CoC_stair6, Deg_stair6, facade_stair6_R, levelSaxis_f, levelPstart(4), stairColu_num, stairL, stairW, CAR, OFFICE, ROOF, Arc_itvl);
[iNO, iEL] = MGT_stair_platform(fileID, iNO, iEL, CoC_stair6, Deg_stair6, levelSaxis, stairColu_num, stairL, stairW, stairB, CAR, OFFICE, ROOF);

fprintf(fileID,'; 塔7\n');
[iNO, iEL] = MGT_side(fileID, iNO, iEL, CoC_side7, facade_side7_R, levelSDaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF, 7, Arc_itvl);

fprintf(fileID,'; 塔8\n');
[iNO, iEL] = MGT_side8(fileID, iNO, iEL, CoC_side8, facade_side8_R, levelSDaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF, Arc_itvl);
fprintf(fileID,'; 塔9\n');
[iNO, iEL] = MGT_side(fileID, iNO, iEL, CoC_side9, facade_side9_R, levelSDaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF, 9, Arc_itvl);
% [iNO, iEL] = MGT_side9(fileID, iNO, iEL, CoC_side9, facade_side9_R, levelSDaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF);

fprintf(fileID,'; 塔10\n');
[iNO, iEL] = MGT_side(fileID, iNO, iEL, CoC_side10, facade_side10_R, levelSDaxis, levelPstart(5), Roof_boundary, CAR, OFFICE, ROOF, 10, Arc_itvl);
%%


%%
fprintf(fileID,'\n');

end