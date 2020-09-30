function r = closestpoint(x, s)
% CLOSESTPOINT - return closest points on line segment s to point x
%   assumes s is a line segment defined by two points
    u = [(s(3)-s(1)) (s(4)-s(2))];
    dotOne = dot(x - [s(1) s(2)], u);   
    dotTwo = dot(u,u);
    t = clamp(dotOne./dotTwo);
    r = (1-t)*([s(1) s(2)]) + t*([s(3) s(4)]);
end

% performs c-like clamp function
function y = clamp(v)
    y = min(max(v, 0), 1);
end
