%% function
% length from P0 to perpendicular point
%
% Xu Yi, 2018

%%
function  L_P = coorPerpL(P0, P1, P2)    % P0��ĳ��,P1\P2�Ǹ�������
% �������㣬���ĳ�㵽���ڸ��������ߵĴ����ľ���
P_m = coorPerp(P0, P1, P2);
L_P = sqrt( (P0(1)-P_m(1))^2 + (P0(2)-P_m(2))^2);
end