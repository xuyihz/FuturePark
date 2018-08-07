%% function
% length from P0 to perpendicular point
%
% Xu Yi, 2018

%%
function  L_P = coorPerpL(C0, P1, P2)    % P0是某点,P1\P2是给定两点
% 给定两点，求出某点到关于该两点连线的垂足点的距离
P_m = coorPerp(C0, P1, P2);
L_P = sqrt( (C0(1)-P_m(1))^2 + (C0(2)-P_m(2))^2);
end