function [y1] = closestpoints(xv, s) % s = [ x1 y1 x2 y2], s[0] = [s(1) s(2)] s[1] = [s(3) s(4)]
% CLOSESTPOINTS - return closest points on line segment s to points xv
%   
    u = [(s(3)-s(1)) (s(4)-s(2))];
    
    ulong = repmat(u, size(xv,1), 1);
    slong = repmat([s(1) s(2)], size(xv,1), 1);
    
    dotOne = dot(xv - slong, ulong, 2);             % dot(xv - slong, ulong);   
    dotTwo = dot(ulong, ulong, 2);                  % f
    
    t = clampVector(dotOne./dotTwo, 0, 1);
    y1 = (1-t)*([s(1) s(2)]) + t*([s(3) s(4)]);
end

% performs c-like clamp function
function y = clampVector(v, lo, hi)
    y = v;
    y(y>hi) = hi;
    y(y<lo) = lo;
end

