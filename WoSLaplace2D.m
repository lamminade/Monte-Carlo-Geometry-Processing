%%%%%%%%%%%%%%%% data %%%%%%%%%%%%%%%%

% a segment is a pair of points
% the scene is an array(vector) of segments
scene = [
    [ [0.5, 0.1], [0.9, 0.5] ] ;
    [ [0.5, 0.9], [0.1, 0.5] ] ;
    [ [0.1, 0.5], [0.5, 0.1] ] ;
    [ [0.5, 0.33333333], [0.5, 0.6666666] ] ;
    [ [0.33333333, 0.5], [0.6666666, 0.5] ] ;
];
   
%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%
out = zeros(129, 129);
s = 129;    % image size
for j = 1:s
    fprintf("row %i of %i\n", j, s);
    for i = 1:s
        x0 = [ (i-1)./(s-1), (j-1)./(s-1) ];
        u = solve(x0, scene);
        out(i,j) = cast(u, 'double');
    end
end

imagesc(out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% returns the point on segment s closest to x
function y1 = closestPoint(x, s) % s = [ x1 y1 x2 y2], s[0] = [s(1) s(2)] s[1] = [s(3) s(4)]
    u = [(s(3)-s(1)) (s(4)-s(2))];
    dotOne = dot(x - [s(1) s(2)], u);   
    dotTwo = dot(u,u);
    t = clamp(dotOne./dotTwo);
    y1 = (1-t)*([s(1) s(2)]) + t*([s(3) s(4)]);
end

% random number in a range - RAND_MAX is like the C variable
function r = myRandom(min, max)
    RAND_MAX = 32767.0;
    rRandMax = 1.0 ./ RAND_MAX;
    u = rRandMax .* (rand() .* RAND_MAX);
    r = (u .* (max - min)) + min;
end

% solves a Laplace equation Δu = 0 at x0, where the boundary is given
% by a collection of segments, and the boundary conditions are given
% by a function g that can be evaluated at any point in space
function res = solve(x0, segments)
    eps = 0.01;     % stopping tolerance
    nWalks = 128;   % number of Monte Carlo samples
    maxSteps = 16;  % maximum walk length
    
    sum = 0;
    for i = 0:nWalks
        x = x0;
        steps = 0;
        
        while true % this is my do-while in matlab
            R = realmax();
            for j = 1:size(segments, 1)                     % # of rows
                % we want the row j of segments
                p = closestPoint(x, segments(j:j,:));   
                R = min(R, norm((x-p)));
            end
            theta = myRandom(0, 2.*pi);
            x = x + [R.*cos(theta), R.*sin(theta)];
            steps = steps + 1;
            
            if (R < eps || steps > maxSteps)
                break
            end
        end
        
        sum = sum + checker(x);
    end
    res = sum ./ nWalks; % monte carlo estimate
end

function c = checker(x)
    s = 12;
    c = fmod(floor(s .* x(1)) + floor(s .* x(2)), 2 );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 helpers for porting                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% performs c-like clamp function
function y = clamp(v)
    y = min(max(v, 0), 1);
end

% performs c-like fmod
function m = fmod(a, b)
    % Where the mod function returns a value in region [0, b), this
    % function returns a value in the region [-b, b), with a negative
    % value only when a is negative.
    
    if a == 0
        m = 0;
    else
        m = mod(a, b) + (b .* (sign(a) - 1) ./ 2);
    end
end

