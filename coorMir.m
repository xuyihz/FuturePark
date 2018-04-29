%% function
% Mirror point
%
% Xu Yi, 23rd April 2018

%%
function  P0_M = coorMir(P0, P1, P2)    % P0是某点,P1\P2是给定两点
% 给定两点，求出某点关于该两点连线的镜像点
if P1(1) == P2(1)   % 关于平行于Y轴的线镜像
    P0_M(1) = P1(1)*2 - P0(1);
    P0_M(2) = P0(2);
elseif P1(2) == P2(2)   % 关于平行于X轴的线镜像
    P0_M(1) = P0(1);
    P0_M(2) = P1(2)*2 - P0(2);
else % 任意点
    % 垂足点
    P_m = coorPerp(P0, P1, P2);
    % 镜像点
    P0_M(1) = P_m(1)*2 - P0(1);
    P0_M(2) = P_m(2)*2 - P0(2);
end
end