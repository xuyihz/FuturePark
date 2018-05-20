%% function
% Coordinate transformation
%
% Xu Yi, 25th March 2018

%%
function  XY_trans = coorTrans(XY, theta)
% 坐标系转换 theta为顺时针为正
% x' = xcos(a) + ysin(a);
% y' = ycos(a) - xsin(a);
XY_trans(1) = XY(1)*cos(theta) + XY(2)*sin(theta);
XY_trans(2) = XY(2)*cos(theta) - XY(1)*sin(theta);
end