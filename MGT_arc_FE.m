%% function
% MGT generate Finite Element(line) to approximate arc
% P_start. P_end, counter-clockwise
% Xu Yi, 2018

%%
function [iNO_end, Deg_num] = MGT_arc_FE(fileID, iNO_start, levelZaxis, CoC, P_start, P_end, interval)
% iNO_end:�սڵ�ű��ݡ�Deg_num:�Ƕȷָ���������ֱ�߶�����
% fileID:MGT�ļ�ID��iNO_start:ʼ�ڵ�š�levelZaxis:�ڵ�Z����ֵ���˴�ΪXoYƽ�满�ߡ�
% CoC:Բ�����ꡣP_start:Բ��������ꡣP_end:Բ���յ����ꡣinterval:�������ֱ�߽��Ƶ�Բ���εĳ��ȡ�
R = sqrt( ( P_start(1)-CoC(1) )^2 + ( P_start(2)-CoC(2) )^2 );  % �뾶
Deg = coorDeg(CoC, P_start, P_end); % Բ���ĽǶ�

Deg_num = 1; Deg_i = Deg/Deg_num;    % ��ʼ����Deg_num:�Ƕȷָ�����Deg_i:������Ԫ�ĽǶ�
while Deg_i*R > interval % ��Deg_i*R > interval��������ָ�
    Deg_num = Deg_num+1;
    Deg_i = Deg/Deg_num;
end

iNO = iNO_start;
if Deg_num == 1 % ˵����Բ��û�зָ������
else
    for i = 1:(Deg_num-1) % Deg_num-1Ϊ�����Ľڵ���
        iNO = iNO+1;
        Deg_temp = -Deg_i*i; % �м�ڵ�����������ת�ĽǶȣ�����coorTransLoc��˳ʱ��Ϊ����coorDeg����ʱ��Ϊ�����������и����š�
        XYcoor = coorTransLoc(CoC, P_start, Deg_temp); % thetaΪ˳ʱ��Ϊ������ֹ�ض�����ۻ�������P_start��ʼ��
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% �ڵ��Ź�����ʱ�롣
            iNO,XYcoor(1),XYcoor(2),levelZaxis);   %
    end
end
iNO_end = iNO;
end
