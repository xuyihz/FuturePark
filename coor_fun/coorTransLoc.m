%% function
% on site Coordinate transformation
%
% Xu Yi, 2018

%%
function  P2 = coorTransLoc(C0, P1, theta) % ԭ������ϵת�� thetaΪ˳ʱ��Ϊ��
% C0Ϊ�ǵ� P1Ϊ��ת���� thetaΪP1��C0ת���ĽǶ� P2ΪP1ת����ĵ�
temp = P1 - C0;
P2 = coorTrans(temp, theta);
P2 = P2 + C0;
end