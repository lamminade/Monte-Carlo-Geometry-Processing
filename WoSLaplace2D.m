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
N = 129;    % image size
out = zeros(N,N);
out2 = zeros(N,N);

for j = 1:N
    fprintf("row %i of %i\n", j, N);
    for i = 1:N
        x0 = [ (i-1)./(N-1), (j-1)./(N-1) ];
        u = laplaceSolve(x0, scene);
        out(i,j) = u;
    end
end

imagesc(out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y1 = closestPoints(xv, s) % s = [ x1 y1 x2 y2], s[0] = [s(1) s(2)] s[1] = [s(3) s(4)]
% closestPoints - 
% returns an array of points on s which are closest
% to the values on x -- actually probably didnt need to vectorize this
% whoops...
%
    u = [(s(3)-s(1)) (s(4)-s(2))];
    
    ulong = repmat(u, size(xv,1), 1);
    slong = repmat([s(1) s(2)], size(xv,1), 1);
    
    dotOne = dot(xv - slong, ulong, 2);             % dot(xv - slong, ulong);   
    dotTwo = dot(ulong, ulong, 2);                  % f
    
    t = clampVector(dotOne./dotTwo, 0, 1);
    y1 = (1-t)*([s(1) s(2)]) + t*([s(3) s(4)]);
end


function res = laplaceSolve(x0, segments)
% laplaceSolve  Laplace Method
% Solves a Laplace equation Δu = 0 at x0, where the boundary is given
%   by a collection of segments, and the boundary conditions are 
%   evaluated by the function 'checker' which can be evaluated at any point in space
%
    eps = 0.01;     % stopping tolerance
    nWalks = 128;   % number of Monte Carlo samples
    maxSteps = 16;  % maximum walk length
    
    % can we do all walks at once...   need to vectorize some of the helper
    % functions.
   
    xv = ones(nWalks,1) * x0;
    steps = 0;
    while true
        R = ones(nWalks,1) * realmax();
        for j = 1:size(segments,1)
            % get closest points simulateously!
            pv = closestPoints(xv, segments(j:j,:));
            R = min(R, vecnorm((xv-pv).').');
        end
        theta = myRandomVect(nWalks, 0, 2.*pi);
        xv = xv + [R.*cos(theta), R.*sin(theta)];
        steps = steps + 1;
        
        if (steps > maxSteps)
            break
        end
    end
    
    sumMat = arrayfun(@(x,y) checker([x y]), xv(:,1:1), xv(:,2:2)); %% or other boundary condition thing
    res = sum(sumMat) ./ nWalks; % monte carlo estimate
end

% boundary condition checker
% the one from the sample code 
function c = checker(x)
    s = 6;
    c = fmod(floor(s .* x(1)) + floor(s .* x(2)), 2 );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 helpers for porting                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% performs c-like clamp function
function y = clampVector(v, lo, hi)
    y = v;
    y(y>hi) = hi;
    y(y<lo) = lo;
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

function rv = myRandomVect(n, min, max)
    RAND_MAX = 32767;
    rRandMax = 1 ./ RAND_MAX;
    u = rRandMax * (rand(n, 1) * RAND_MAX);
    rv = (u * (max - min)) + min;
end

