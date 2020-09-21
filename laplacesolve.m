function res = laplacesolve(x0, segments, boundaryConditionFxn)
% LAPLACESOLVE - Solves a Laplace equation Δu = 0 at x0
% Solves a Laplace equation Δu = 0 at x0, where the boundary is given
%   by a collection of segments, and the boundary conditions are 
%   evaluated by the function 'checker' which can be evaluated at any point in space

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
            pv = closestpoints(xv, segments(j:j,:));
            R = min(R, vecnorm((xv-pv).').');
        end
        theta = randomvector(nWalks, 0, 2.*pi);
        xv = xv + [R.*cos(theta), R.*sin(theta)];
        steps = steps + 1;
        
        if (steps > maxSteps)
            break
        end
    end
    
    sumMat = arrayfun(@(x,y) boundaryConditionFxn([x y]), xv(:,1:1), xv(:,2:2)); %% or other boundary condition thing
    res = sum(sumMat) ./ nWalks; % monte carlo estimate
end

function r = randomvector(n, min, max)
% RANDOMVECTOR - returns a (n 1) vector of random numbers in a range
%   returns a vector matrix of dimension n x 1 using the rand function to 
% draw the values from a uniform distribution in the open interval, (min, max).
  r = (max - min).*rand(n,1) + min;
end