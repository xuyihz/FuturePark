%% function
% MGT initial conditions
% 
% Xu Yi, 19th March 2018
% Xu Yi, 23rd March 2018, revised

%%
function MGT_init(fileID)
%%
fprintf(fileID,...
    ';---------------------------------------------------------------------------\n'...
    );
fprintf(fileID,';  midas Gen Text(MGT) File.\n');
fprintf(fileID,';  Date : %s\n',datetime('today'));
fprintf(fileID,...
    ';---------------------------------------------------------------------------\n'...
    );
fprintf(fileID,'\n');

%% VERSION
fprintf(fileID,'*VERSION\n');
fprintf(fileID,'   8.6.5\n');
fprintf(fileID,'\n');

%% UNIT
FORCE = 'KN'; LENGTH = 'MM'; HEAT = 'KJ'; TEMPER = 'C';

fprintf(fileID,'*UNIT    ; Unit System\n');
fprintf(fileID,'; FORCE, LENGTH, HEAT, TEMPER\n');
fprintf(fileID,'   %s, %s, %s, %s\n',FORCE,LENGTH,HEAT,TEMPER);
fprintf(fileID,'\n');

%% STRUCTYPE
iSTYP = '0'; iMASS = '1'; iSMAS = '1'; bMASSOFFSET = 'NO';
bSELFWEIGHT = 'YES'; GRAV = '9806'; TEMPER = '0';
bALIGNBEAM = 'NO'; bALIGNSLAB = 'NO'; bROTRIGID = 'NO';

fprintf(fileID,'*STRUCTYPE    ; Structure Type\n');
fprintf(fileID,'; iSTYP, iMASS, iSMAS, bMASSOFFSET, bSELFWEIGHT, GRAV, TEMPER, bALIGNBEAM, bALIGNSLAB, bROTRIGID\n');
fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
    iSTYP, iMASS, iSMAS, bMASSOFFSET, bSELFWEIGHT, ...
    GRAV, TEMPER, bALIGNBEAM, bALIGNSLAB, bROTRIGID);
fprintf(fileID,'\n');

%% REBAR-MATL-CODE
CONC_CODE = 'GB10(RC)'; CONC_MDB = 'HRB400';
SRC_CODE = 'GB10(RC)'; SRC_MDB = 'HRB400';

fprintf(fileID,'*REBAR-MATL-CODE    ; Rebar Material Code\n');
fprintf(fileID,'; CONC_CODE, CONC_MDB, SRC_CODE, SRC_MDB\n');
fprintf(fileID,'   %s, %s, %s, %s\n',CONC_CODE,CONC_MDB,SRC_CODE,SRC_MDB);
fprintf(fileID,'\n');

%% STLDCASE
LCNAME = {'DL';'LL';'WX';'WY'};
LCTYPE = {'D';'L';'W';'W'};
DESC = {'Dead Load';'Live Load';'Wind X';'Wind Y'};

fprintf(fileID,'*STLDCASE    ; Static Load Cases\n');
fprintf(fileID,'; LCNAME, LCTYPE, DESC\n');
for i = 1:length(LCNAME)
    LCN = LCNAME{i}; LCT = LCTYPE{i}; DES = DESC{i};
    fprintf(fileID,'   %s, %s, %s\n',LCN,LCT,DES);
end
fprintf(fileID,'\n');

end