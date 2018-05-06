function [NURBS_curve, N, Rational_func, S, U] = nurbsfun(NURBS_order, knotV, weightV, control_pts, U)
% Evaluate NURBS at specified locations.
%
% Input arguments:
% NURBS_order:
%    NURBS order (2 for linear, 3 for quadratic, 4 for cubic, etc.)
% knotV:
%    knot vector
% weightV:
%    weight vector
% control_pts:
%    control points, typically 2-by-m
% U (optional):
%    values where the NURBS is to be evaluated, or a positive
%    integer to set the number of points to automatically allocate
%
% Output arguments:
% NURBS_curve:
%    points of the NURBS curve
% N:
%    basis function value for each element in U
% Rational_func:
%    rational function value for each element in U
% S:
%    column containing the value of the summarize of basis function per
%    weight
% U:
%    points where the NURBS curve is evaluated

%Written by Graziano Fuccio, email: g.fuccio359@gmail.com
validateattributes(NURBS_order, {'numeric'}, {'positive','integer','scalar'});
assert(all( knotV(2:end)-knotV(1:end-1) >= 0 ), 'nurbs:InvalidArgumentValue', ...
    'Knot vector values should be nondecreasing.');
validateattributes(control_pts, {'numeric'}, {'real','2d'});
nctrl = numel(knotV)- NURBS_order;
assert(size(control_pts,2) == nctrl, 'nurbs:DimensionMismatch', ...
    'Invalid number of control points, %d given, %d required.', size(control_pts,2), nctrl);
assert(size(control_pts,2) == numel(weightV), 'nurbs:DimensionMismatch', ...
    'Invalid number of weight points, %d given, %d required.', numel(weightV), size(control_pts,2));
if nargin < 5
    U = linspace(knotV(NURBS_order), knotV(end-NURBS_order), 10*size(control_pts,2));% allocate points uniformly,
elseif isscalar(U) && U > 1
    validateattributes(U, {'numeric'}, {'positive','integer','scalar'});
    U = linspace(knotV(NURBS_order), knotV(end-NURBS_order), U);  % allocate points uniformly
end

totalU = numel(U);%total elements to evaluate
nc = size(control_pts, 2);
N = zeros(totalU, nc);

%calcolate basis function
for i = 1 : totalU
    u = U(i);
    N(i, :) = basisfunction(NURBS_order, nc, u, knotV);
end

%calcolate denominator of rational basis functions
S = zeros(totalU, 1);
for i = 1 : totalU
    tmp = N(i, :);
    for j = 1 : size(tmp, 2)
        S(i) = S(i) + (tmp(j) * weightV(j));
    end
end

%calcolate rational basis functions
Rational_func = zeros(totalU, nc);
for i = 1 : totalU
    tmp = N(i, :);
    for j = 1 : size(tmp, 2)
        if(S(i) ~= 0)
            Rational_func(i, j) = (tmp(j) * weightV(j)) / S(i);
        else
            Rational_func(i, j) = 0;
        end
    end
end

%calcolate curve points
NURBS_curve = zeros(2, totalU);
for i = 1 : totalU
    tmp = Rational_func(i, :);
    sum = [0; 0];
    for j = 1 : size(tmp, 2)
        sum = sum + (control_pts(:, j) * tmp(j));
    end
    NURBS_curve(:, i) = sum;
end
end