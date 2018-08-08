%% function
% length from P0 to perpendicular point
%
% Xu Yi, 2018

%%
function  L_P = coorPerpL(C0, P1, P2)    % P0��ĳ��,P1\P2�Ǹ�������
% �������㣬���ĳ�㵽���ڸ��������ߵĴ����ľ���
P_m = coorPerp(C0, P1, P2);
L_P = sqrt( (C0(1)-P_m(1))^2 + (C0(2)-P_m(2))^2);
end