%% function
% MGT generate Finite Element(line) to approximate arc
% P_start. P_end, counter-clockwise
% Xu Yi, 2018

%%
function [iNO_end, Deg_num] = MGT_arc_FE(fileID, iNO_start, levelZaxis, CoC, P_start, P_end, interval)
% iNO_end:终节点号备份。Deg_num:角度分隔总数，即直线段数。
% fileID:MGT文件ID。iNO_start:始节点号。levelZaxis:节点Z坐标值，此处为XoY平面弧线。
% CoC:圆心坐标。P_start:圆弧起点坐标。P_end:圆弧终点坐标。interval:最大容许直线近似的圆弧段的长度。
R = sqrt( ( P_start(1)-CoC(1) )^2 + ( P_start(2)-CoC(2) )^2 );  % 半径
Deg = coorDeg(CoC, P_start, P_end); % 圆弧的角度

Deg_num = 1; Deg_i = Deg/Deg_num;    % 初始化。Deg_num:角度分隔数，Deg_i:单个单元的角度
while Deg_i*R > interval % 如Deg_i*R > interval，则继续分割
    Deg_num = Deg_num+1;
    Deg_i = Deg/Deg_num;
end

iNO = iNO_start;
if Deg_num == 1 % 说明此圆弧没有分割，故跳过
else
    for i = 1:(Deg_num-1) % Deg_num-1为新增的节点数
        iNO = iNO+1;
        Deg_temp = -Deg_i*i; % 中间节点相对于起点旋转的角度，由于coorTransLoc是顺时针为正，coorDeg是逆时针为正，故这里有个负号。
        XYcoor = coorTransLoc(CoC, P_start, Deg_temp); % theta为顺时针为正。防止截断误差累积，都从P_start开始。
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...	% 节点编号规则：逆时针。
            iNO,XYcoor(1),XYcoor(2),levelZaxis);   %
    end
end
iNO_end = iNO;
end
