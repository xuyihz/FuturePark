%control points
P = [5 2 2 6 10 10 7; 1 3 8 10 8 3 1];
n = 3; %curve's degree n = p + 1 where p is the degree of the curve
%knot vector
t = [0 0 0 0.1 0.2 0.5 0.7 1 1 1];
%weight vector
w = [1 1 1 1 1 1 1];
%call of the function
C = nurbsfun(n, t, w, P);
hold on
xlim([0, 12])
ylim([0, 12])
plot(P(1, :), P(2, :), 'bo')
plot(P(1, :), P(2, :), 'g')
plot(C(1, :), C(2, :), 'r')
%changing weight
w = [1 1 1 3 1 2 0];
C = nurbsfun(n, t, w, P);
plot(C(1, :), C(2, :), 'b')
%cambio pesi
w = [1 1 3 1 4 1 1];
C = nurbsfun(n, t, w, P);
plot(C(1, :), C(2, :), 'm')