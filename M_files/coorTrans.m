%% function
% Coordinate transformation
%
% Xu Yi, 25th March 2018

%%
function  [x_trans,y_trans] = coorTrans(x, y, theta)
% ����ϵת�� thetaΪ˳ʱ��Ϊ��
% x' = xcos(a) + ysin(a);
% y' = ycos(a) - xsin(a);
x_trans = x*cos(theta) + y*sin(theta);
y_trans = y*cos(theta) - x*sin(theta);
end