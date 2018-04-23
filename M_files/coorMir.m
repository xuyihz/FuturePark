%% function
% Mirror point
%
% Xu Yi, 23rd April 2018

%%
function  P0_M = coorMir(P0, P1, P2)
% �������㣬���ĳ����ڸ��������ߵľ����
if P1(1) == P2(1)   % ����ƽ����Y����߾���
    P0_M(1) = P1(1)*2 - P0(1);
    P0_M(2) = P0(2);
elseif P1(2) == P2(2)   % ����ƽ����X����߾���
    P0_M(1) = P0(1);
    P0_M(2) = P1(2)*2 - P0(2);
else % �����
    a = P2(2) - P1(2);
    b = P1(1) - P2(1);
    c = P2(1) * P1(2) - P1(1) * P2(2);
    sqr = a * a + b * b;
    % �����
    P_m(1) = (b * b * P0(1) - a * b * P0(2) - a * c) / sqr;
    P_m(2) = (a * a * P0(2) - a * b * P0(1) - b * c) / sqr;
    % �����
    P0_M(1) = P_m(1)*2 - P0(1);
    P0_M(2) = P_m(2)*2 - P0(2);
end
end