%% function
% 3 points, get its degree
%
% Xu Yi, 2018

%%
function  Deg = coorDeg(C0, P1, P2)    % ��֪���㣬 �õ���C0��Ϊ�ǵ�ļн� % C0Բ��/P1P2����
% ��֪���㣬 �õ���C0��Ϊ�ǵ�ļнǡ� P1->P2��ʱ��
% ��������������������ģ�Ļ�����������нǵ�����

C0P1 = P1 - C0; % ����1
C0P2 = P2 - C0; % ����2

C0P1_n = norm(C0P1);
C0P2_n = norm(C0P2);
C0P1P2 = dot(C0P1, C0P2);

Deg = acos( C0P1P2 / ( C0P1_n*C0P2_n ) );
end