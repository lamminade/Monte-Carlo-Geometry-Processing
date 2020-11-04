function res = poissonsolve(x0, segments, f, g)
% POISSONSOLVE - same as laplacesolve except with poisson sampling of
% another point y
% solves a Laplace equation Δu = f at x0, where the boundary is given
% by a collection of segments, and the boundary conditions are given
% by a function g that can be evaluated at any point in space

    eps = 0.001;   % stopping tolerance
    nWalks = 32;   % num Monte carlo samples
    maxSteps = 16; % maximum walk length
    
    sumv = zeros(nWalks, 1);
    xv = ones(nWalks,1) * x0;
    steps = 0;
    while true
        R = ones(nWalks,1) * realmax();
        for j = 1:size(segments,1)
            % get closest points simulateously!
            pv = closestpoints(xv, segments(j:j,:));
            R = min(R, vecnorm((xv-pv).').');
        end
            
        % sample points y uniformly from the ball of radius R around x
        rVect = R .* sqrt(randomvector(nWalks, 0, 1));
        alpha = randomvector(nWalks, 0, 1);
        yVect = xv + [rVect .* cos(alpha), rVect .* sin(alpha)];
        sumv = sumv + (pi .* R .* R) .* f(yVect) .* harmonicgreensv(rVect, R);
            
        % sample the next point x uniformly from the sphere around x
        theta = randomvector(nWalks, 0, 2 * pi);
        xv = xv + [R.*cos(theta), R.*sin(theta)];
        steps = steps + 1;

        if (steps > maxSteps)
            break;
        end
    end
    
    % apply g to all rows of xv
    gxv = arrayfun(@(x,y) g([x y]), xv(:,1:1), xv(:,2:2));
    sumv = sumv + gxv;
    res = sum(sumv) ./ nWalks; % Monte Carlo estimate
end

function GrR = harmonicgreens(r, R)
% harmonic Green's function 
% for a 2D ball of radius R
	GrR = log(R./r) / (2 * pi);
    GrR(GrR == -inf) = 0;
end

function r = randomvector(n, min, max)
% RANDOMVECTOR - returns a (n 1) vector of random numbers in a range
%   returns a vector matrix of dimension n x 1 using the rand function to 
% draw the values from a uniform distribution in the open interval, (min, max).
  r = (max - min).*rand(n,1) + min;
end