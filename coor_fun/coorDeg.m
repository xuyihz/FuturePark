%% function
% 3 points, get its degree
%
% Xu Yi, 2018

%%
function  Deg = coorDeg(C0, P1, P2)    % 已知三点， 得到以C0点为角点的夹角 % C0圆心/P1P2两点
% 已知三点， 得到以C0点为角点的夹角。 P1->P2逆时针
% 向量的数量积除以向量模的积等于向量间夹角的余弦

C0P1 = P1 - C0; % 向量1
C0P2 = P2 - C0; % 向量2

C0P1_n = norm(C0P1);
C0P2_n = norm(C0P2);
C0P1P2 = dot(C0P1, C0P2);

Deg = acos( C0P1P2 / ( C0P1_n*C0P2_n ) );
end