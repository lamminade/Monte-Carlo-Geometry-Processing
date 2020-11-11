scene= [
    [ [50, 50], [50, 50] ] ;
    [ [0, 0], [1, 0] ] ;
    [ [0, 0], [0, 1] ] ;
    [ [0, 1], [1, 1] ] ;
    [ [1, 0], [1, 1] ] ;
];          % need some scene
c = 1;      % need some c


% To validate the implementation we solve the Poisson equation
%    
%   Δu - cu = Δu0  on Ω
%	      u =  u0  on ∂Ω
%    
% where u0 is some reference function and c is some constant.
% The solution should converge to u = u0 as the number of samples
% N increases.
N = 100; % 28;
out = zeros(N,N);
for j = 1:N
	for i = 1:N
        x0 = [ i/N j/N ];
        % solves a Laplace equation Δu = f at x0, where the boundary is given
        % by a collection of segments, and the boundary conditions are given
        % by a function g that can be evaluated at any point in space
        % screenedpoissonsolve(x0, segments, f, g, c)
        u = screenedpoissonsolve(x0, scene, @laplace_urefv, @uref, c);
        out(i,j) = u;
	end
end

disp("fin");
imagesc(out);

% reference solution
function c = uref(x)
   c = cos(2 * pi * x(1)) * sin(2 * pi * x(2));
end

% Laplacian of reference solution
function cv = laplace_urefv(xv)
   cv = 8 .* (pi .* pi) .* cos(2.*pi.*xv(:,1:1)) .* sin(2 .* pi .* xv(:,2:2));
end