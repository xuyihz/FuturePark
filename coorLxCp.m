%% function
% line x circle point
%
% Xu Yi, 26th April 2018

%%
function  [P_LxC1,P_LxC2] = coorLxCp(P0, P1, P2, R)    % 线与圆交点 % 仅限特定情况。P0=[0,0]
% 给定圆心P0半径R，和另两点P1\P2，求出两点连线与圆交点的其中一点。(未考虑特殊情况，即两点连线平行与坐标轴情况。)
P_m = coorPerp(P0, P1, P2); % 得到垂足点
P_LxC(1) = -sqrt(P_m(1)^2+P_m(2)^2);	% 假设直线平行于Y轴情况，P_LxC(1)应等于-P_m到P0距离(特殊情况)
P_LxC(2) = -sqrt(R^2-P_LxC(1)^2);   % P_LxC(1)^2+P_LxC(2)^2=R^2
theta = atan( P_m(2)/P_m(1) );
[P_LxC1,P_LxC2] = coorTrans(P_LxC(1), P_LxC(2), -theta);    % 逆时针转，故有负号
end