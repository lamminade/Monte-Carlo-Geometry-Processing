scene= [
    [ [50, 50], [50, 50] 1 ] ;
    [ [0, 0], [1, 0] -1 ] ;
    [ [0, 0], [0, 1]  -1 ] ;
    [ [0, 1], [1, 1]  -1 ] ;
    [ [1, 0], [1, 1]  -1 ] ;
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

% reference solution (g)
function c = uref(x)
	c = cos(2 * pi * x(1)) * sin(2 * pi * x(2)); 
%     if (x(2) >= 40)
%         c = 10;
%     end
%     if (x(2) < 40)
%         c = 1;
%     end
end

% Laplacian of reference solution (f)
function c = laplace_urefv(xv)
    c = 8 .* (pi .* pi) .* cos(2.*pi.*xv(:,1:1)) .* sin(2 .* pi .* xv(:,2:2));
%     [h, w] = size(segments);
%     
%     distances = zeros(size(xv,1), h); 
%     for j = 1:h
%     	pv = closestpoints(xv, segments(j:j,:));
%         distances(:,j:j) = vecnorm((xv-pv).').';
%     end
%     
%     [~, I] = min(distances,[],2);
%     
%     c = zeros(size(I,1),1);
%     for j = 1:size(I,1)
%         % how do we want to treat the value of the boundary 
%         s = segments(I(j),w);
%         % just set it = to boundary condition?
%         c(j,1) = s; 
%     end
%     if (xv(2) >= 40)
%         cv = 10;
%     end
%     if (xv(2) < 40)
%         cv = 1;
%     end
    %xv(xv(:,2) > 40) = 10;
    %xv(xv(:,2) < 40) = 1;
end