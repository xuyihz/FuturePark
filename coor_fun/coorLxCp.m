%% function
% line x circle point
%
% Xu Yi, 2018

%%
function  [X, Y, Len] = coorLxCp(C0, R, P1, P2)    % ����Բ����
% ����Բ��P0�뾶R����������P1\P2���������������Բ��������㡣
L_P = coorPerpL(C0, P1, P2); % ���߳�
if L_P > R % ֱ����Բ���ཻ
    X = 0; Y = 0; Len = 0;
else
    P_m = coorPerp(C0, P1, P2); % �õ������
    if L_P == R % ֱ����Բ����
        X = P_m(1); Y = P_m(2); Len = 1;
    else % ֱ����Բ�ཻ
        Line_L = sqrt(R^2 - L_P^2); % �����������һ�ߵı߳�
        if L_P == 0 % ֱ�߾���Բ��
            if P1(1) == P2(1)   % ֱ��ƽ����Y��
                X(1) = C0(1); X(2) = C0(1);
                Y(1) = C0(2) + R; Y(2) = C0(2) - R;
            elseif P1(2) == P2(2)   % ֱ��ƽ����X��
                X(1) = C0(1) + R; X(2) = C0(1) - R;
                Y(1) = C0(2); Y(2) = C0(2);
            else % һ��ֱ��
                P1_P2 = P1 - P2;
                P1_P2_theta = arctan( P1_P2(2)/P1_P2(1) );
                X_delta = R * cos(P1_P2_theta); % 
                Y_delta = R * sin(P1_P2_theta);
                if P1_P2(1)*P1_P2(2) > 0 % ֱ��б��Ϊ��
                    X(1) = C0(1) + X_delta; X(2) = C0(1) - X_delta;
                    Y(1) = C0(2) + Y_delta; Y(2) = C0(2) - Y_delta;
                else
                    X(1) = C0(1) + X_delta; X(2) = C0(1) - X_delta;
                    Y(1) = C0(2) + Y_delta; Y(2) = C0(2) - Y_delta;
                end
            end
        else % ֱ�߲�����Բ��
            X_delta = Line_L / L_P * ( P_m(2)-C0(2) ); % ���Ƕ�Ӧ��ϵ
            Y_delta = Line_L / L_P * ( P_m(1)-C0(1) );
            X(1) = P_m(1) + X_delta; X(2) = P_m(1) - X_delta;
            Y(1) = P_m(2) - Y_delta; Y(2) = P_m(2) + Y_delta;
        end
        Len = 2; % ��ĸ���
    end
end
end