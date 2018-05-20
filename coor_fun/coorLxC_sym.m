%% function
% line x circle
%
% Xu Yi, 2018

%%
function  [X, Y, Len] = coorLxC_sym(C0, R, P1, P2)    % 线与圆交点 % C0圆心/R半径/P1P2直线上两点
% 给定圆心C0半径R，和另两点P1\P2，求出两点连线与圆交点。
syms x y

eqn1 = (x-C0(1))^2 + (y-C0(2))^2 == R^2;
if P1(1) == P2(1)
    eqn2 = x == P1(1); % 平行于Y轴
elseif P1(2) == P2(2)
    eqn2 = y == P1(2); % 平行于X轴
else
    eqn2 = (y-P1(2))/(P2(2)-P1(2)) == (x-P1(1))/(P2(1)-P1(1));
end

sol = solve( eqn1, eqn2, x,y);

X = double(sol.x);
Y = double(sol.y);
Len = length (X); % 解的个数
end