%% function
% perpendicular point
%
% Xu Yi, 26th April 2018

%%
function  P_m = coorPerp(P0, P1, P2)    % P0��ĳ��,P1\P2�Ǹ�������
% �������㣬���ĳ����ڸ��������ߵĴ����
if P1(1) == P2(1)   % ����ƽ����Y����߾���
    P_m(1) = P1(1);
    P_m(2) = P0(2);
elseif P1(2) == P2(2)   % ����ƽ����X����߾���
    P_m(1) = P0(1);
    P_m(2) = P1(2);
else % ����� % ֱ�߷��� ax+by+c=0
    a = P2(2) - P1(2);
    b = P1(1) - P2(1);
    c = P2(1) * P1(2) - P1(1) * P2(2);
    sqr = a * a + b * b;
    % �����
    P_m(1) = (b * b * P0(1) - a * b * P0(2) - a * c) / sqr;
    P_m(2) = (a * a * P0(2) - a * b * P0(1) - b * c) / sqr;
end
end