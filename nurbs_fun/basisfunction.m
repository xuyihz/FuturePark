function [N] = basisfunction(NURBS_order, npts, u, knotV)
% Evaluate basis function at specified u.
%
% Input arguments:
% NURBS_order:
%    NURBS order (2 for linear, 3 for quadratic, 4 for cubic, etc.)
% npts:
%    number of control points
% u:
%    point to evaluate
% knotV:
%    knot vector
%
% Output arguments:
% N:
%   vector with size npts containing value of the basis function at u

%Written by Graziano Fuccio, email: g.fuccio359@gmail.com
nplusc = npts + NURBS_order;
N = zeros(npts);

for i = 1: nplusc - 1
    if u >= knotV(i) && u <= knotV(i+1)
        N(i) = 1;
    else
        N(i) = 0;
    end
end

%application of the formula that you can find on the NURBS book
for k = 2 : NURBS_order
    for i = 1 : nplusc - k
        if N(i) ~= 0
            d = ((u - knotV(i)) * N(i)) / (knotV(i + k - 1) - knotV(i));
        else
            d = 0;
        end
        if(N(i + 1) ~= 0)
            e = ((knotV(i + k) - u) * N(i + 1)) / (knotV(i + k) - knotV(i + 1));
        else
            e = 0;
        end
        N(i) = d + e;
    end
end

N = N(:, 1);
N = N';

end