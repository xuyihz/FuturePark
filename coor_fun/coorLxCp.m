%% function
% line x circle point
%
% Xu Yi, 26th April 2018

%%
function  [P_LxC1,P_LxC2] = coorLxCp(P0, P1, P2, R)    % ����Բ���� % �����ض������P0=[0,0]
% ����Բ��P0�뾶R����������P1\P2���������������Բ���������һ�㡣(δ�����������������������ƽ���������������)
P_m = coorPerp(P0, P1, P2); % �õ������
P_LxC(1) = -sqrt(P_m(1)^2+P_m(2)^2);	% ����ֱ��ƽ����Y�������P_LxC(1)Ӧ����-P_m��P0����(�������)
P_LxC(2) = -sqrt(R^2-P_LxC(1)^2);   % P_LxC(1)^2+P_LxC(2)^2=R^2
theta = atan( P_m(2)/P_m(1) );
[P_LxC1,P_LxC2] = coorTrans(P_LxC(1), P_LxC(2), -theta);    % ��ʱ��ת�����и���
end