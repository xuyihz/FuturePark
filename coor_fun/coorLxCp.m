%% function
% line x circle point
%
% Xu Yi, 2018

%%
function  [X, Y, Len] = coorLxCp(C0, R, P1, P2)    % 线与圆交点
% 给定圆心P0半径R，和另两点P1\P2，求出两点连线与圆交点的两点。
L_P = coorPerpL(C0, P1, P2); % 垂线长
if L_P > R % 直线与圆不相交
    X = 0; Y = 0; Len = 0;
else
    P_m = coorPerp(C0, P1, P2); % 得到垂足点
    if L_P == R % 直线与圆相切
        X = P_m(1); Y = P_m(2); Len = 1;
    else % 直线与圆相交
        Line_L = sqrt(R^2 - L_P^2); % 求得三角形另一边的边长
        if L_P == 0 % 直线经过圆心
            if P1(1) == P2(1)   % 直线平行于Y轴
                X(1) = C0(1); X(2) = C0(1);
                Y(1) = C0(2) + R; Y(2) = C0(2) - R;
            elseif P1(2) == P2(2)   % 直线平行于X轴
                X(1) = C0(1) + R; X(2) = C0(1) - R;
                Y(1) = C0(2); Y(2) = C0(2);
            else % 一般直线
                P1_P2 = P1 - P2;
                P1_P2_theta = arctan( P1_P2(2)/P1_P2(1) );
                X_delta = R * cos(P1_P2_theta); % 
                Y_delta = R * sin(P1_P2_theta);
                if P1_P2(1)*P1_P2(2) > 0 % 直线斜率为正
                    X(1) = C0(1) + X_delta; X(2) = C0(1) - X_delta;
                    Y(1) = C0(2) + Y_delta; Y(2) = C0(2) - Y_delta;
                else
                    X(1) = C0(1) + X_delta; X(2) = C0(1) - X_delta;
                    Y(1) = C0(2) + Y_delta; Y(2) = C0(2) - Y_delta;
                end
            end
        else % 直线不经过圆心
            X_delta = Line_L / L_P * ( P_m(2)-C0(2) ); % 三角对应关系
            Y_delta = Line_L / L_P * ( P_m(1)-C0(1) );
            X(1) = P_m(1) + X_delta; X(2) = P_m(1) - X_delta;
            Y(1) = P_m(2) - Y_delta; Y(2) = P_m(2) + Y_delta;
        end
        Len = 2; % 解的个数
    end
end
end