%% function
% on site Coordinate transformation
%
% Xu Yi, 2018

%%
function  P2 = coorTransLoc(C0, P1, theta) % 原点坐标系转换 theta为顺时针为正
% C0为角点 P1为待转换点 theta为P1绕C0转换的角度 P2为P1转换后的点
temp = P1 - C0;
P2 = coorTrans(temp, theta);
P2 = P2 + C0;
end