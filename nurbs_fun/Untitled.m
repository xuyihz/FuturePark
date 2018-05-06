%control points
P = [0, 3603, 10809, 3603, 0; 21800, 18167, 10900, 3633, 0];
n = 4; %curve's degree n = p + 1 where p is the degree of the curve
%knot vector
t = [0 0 0 0 0.5 1 1 1 1];
%weight vector
w = ones(1,5);
%call of the function
C = nurbsfun(n, t, w, P);
hold on
plot(P(1, :), P(2, :), 'bo')
plot(P(1, :), P(2, :), 'g')
plot(C(1, :), C(2, :), 'r')