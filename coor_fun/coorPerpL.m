%% function
% length from P0 to perpendicular point
%
% Xu Yi, 2018

%%
function  L_P = coorPerpL(P0, P1, P2)    % P0是某点,P1\P2是给定两点
% 给定两点，求出某点到关于该两点连线的垂足点的距离
P_m = coorPerp(P0, P1, P2);
L_P = sqrt( (P0(1)-P_m(1))^2 + (P0(2)-P_m(2))^2);
end