%% function
% MGT material
% 
% Xu Yi, 19th March 2018

%%
function MGT_mat(fileID)
%%
fprintf(fileID,'*MATERIAL    ; Material\n');
fprintf(fileID,'; iMAT, TYPE, MNAME, SPHEAT, HEATCO, PLAST, TUNIT, bMASS, DAMPRATIO, [DATA1]           ; STEEL, CONC, USER\n; iMAT, TYPE, MNAME, SPHEAT, HEATCO, PLAST, TUNIT, bMASS, DAMPRATIO, [DATA2], [DATA2]  ; SRC\n; [DATA1] : 1, DB, NAME, CODE, USEELAST, ELAST\n; [DATA1] : 2, ELAST, POISN, THERMAL, DEN, MASS\n; [DATA1] : 3, Ex, Ey, Ez, Tx, Ty, Tz, Sxy, Sxz, Syz, Pxy, Pxz, Pyz, DEN, MASS         ; Orthotropic\n; [DATA2] : 1, DB, NAME, CODE, USEELAST, ELAST or 2, ELAST, POISN, THERMAL, DEN, MASS\n');

fprintf(fileID,'   1, STEEL, Q345, 0, 0, , C, NO, 0.02, 1, GB 50917-13(S), , Q345, NO, 206000\n');
fprintf(fileID,'   2, CONC, C30, 0, 0, , C, NO, 0.05, 1, GB 50917-13(RC), , C30, NO, 30000\n');

%%
fprintf(fileID,'\n');

end