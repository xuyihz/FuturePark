%% function
% line x circle
%
% Xu Yi, 2018

%%
function  [X, Y, Len] = coorLxC_sym(C0, R, P1, P2)    % ����Բ���� % C0Բ��/R�뾶/P1P2ֱ��������
% ����Բ��C0�뾶R����������P1\P2���������������Բ���㡣
syms x y

eqn1 = (x-C0(1))^2 + (y-C0(2))^2 == R^2;
if P1(1) == P2(1)
    eqn2 = x == P1(1); % ƽ����Y��
elseif P1(2) == P2(2)
    eqn2 = y == P1(2); % ƽ����X��
else
    eqn2 = (y-P1(2))/(P2(2)-P1(2)) == (x-P1(1))/(P2(1)-P1(1));
end

sol = solve( eqn1, eqn2, x,y);

X = double(sol.x);
Y = double(sol.y);
Len = length (X); % ��ĸ���
end