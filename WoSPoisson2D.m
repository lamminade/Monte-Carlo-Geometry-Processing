% segments enclosing unit square
scene = [
    [ [0, 0], [1, 0] ] ;
    [ [1, 0], [1, 1] ] ;
    [ [1, 1], [0, 1] ] ;
    [ [0, 1], [0, 0] ] ;
];

% To validate the implementation we solve the Poisson equation
%    
%   Δu = Δu0  on Ω
%	u =  u0  on ∂Ω
%    
% where u0 is some reference function. The solution should
% converge to u = u0 as the number of samples N increases.
N = 300;
out = zeros(N,N);
for j = 1:N
	for i = 1:N
        x0 = [ i/N j/N ];
        u = poissonsolve(x0, scene, @laplace_urefv, @uref);
        out(i,j) = u;
	end
end
 
disp("fin");
imagesc(out);

function r = myrandom(min, max)
% MYRANDOM - returns a random value in the range [rMin, rMax]
    RAND_MAX = 32767.0;
    rRandMax = 1.0 ./ RAND_MAX;
    u = rRandMax .* (rand() .* RAND_MAX);
    r = (u .* (max - min)) + min;
end

function GrR = harmonicgreens(r, R)
% HARMONICGREENS - harmonic Green's function for a 2D ball of radius R
	GrR = log(R/r) / (2 * pi);
    if (GrR == -inf)
        GrR = 0;
    end
end

% solves a Laplace equation Δu = f at x0, where the boundary is given
% by a collection of segments, and the boundary conditions are given
% by a function g that can be evaluated at any point in space
function res = poissonsolveOne(x0, segments, f, g)
    eps = 0.001;   % stopping tolerance
    nWalks = 32;   % num Monte carlo samples
    maxSteps = 16; % maximum walk length
    
    sum = 0;
    for i = 0:nWalks
        x = x0;
        R = 0;
        steps = 0;
        while(true) 
            % get the distance to the closest point on any segment
            R = realmax();
            for t = 1:size(size(segments,1))
                p = closestpoint(x, segments(t:t,:));
                R = min(R, norm(x-p));
            end
            
            % sample a point y uniformly from the ball of radius R around x
            r = R * sqrt(myrandom(0,1));            
            alpha = myrandom(0, 1);
            y = x + [r*cos(alpha), r*sin(alpha)];
            sum = sum + (pi * R * R) * f(y) * harmonicgreens(r, R);
            
            % sample the next point x uniformly from the sphere around x
            theta = myrandom(0, 2 * pi);
            x = x + [R*cos(theta), R*sin(theta)];
            steps = steps + 1;
            
            if (R < eps || steps > maxSteps)
                break;
            end
        end
        
        sum = sum + g(x);
    end
    
    res = sum / nWalks; % Monte Carlo estimate
end

% reference solution
function c = uref(x)
   c = cos(2 * pi * x(1)) * sin(2 * pi * x(2));
end

% Laplacian of reference solution
function c = laplace_uref(x)
   c = 8 * (pi * pi) * cos(2 * pi * x(1)) * sin(2 * pi * x(2));
end

% Laplacian of reference solution
function cv = laplace_urefv(xv)
   cv = 8 .* (pi .* pi) .* cos(2.*pi.*xv(:,1:1)) .* sin(2 .* pi .* xv(:,2:2));
end




