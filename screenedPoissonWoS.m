% scene= [
%     [ [50, 50], [50, 50] 1 ] ;
%     [ [0, 0], [1, 0] -1 ] ;
%     [ [0, 0], [0, 1]  -1 ] ;
%     [ [0, 1], [1, 1]  -1 ] ;
%     [ [1, 0], [1, 1]  -1 ] ;
% ];          % need some scene

scene = [
     [ [0.6, 0.49], [0.6, 0.5] 1]; % top wall
%     [ [0.7, 0.45], [0.7, 0.7] 1]; % bottom wall (has opening)
%     [ [0.3, 0.3], [0.7, 0.3] 1]; % left wall
%     [ [0.3, 0.7], [0.7, 0.7] 1]; % right wall 
     % outside perimeter of image ! 
%     [ [0, 0], [0, 1] -1];    % LS
%     [ [1, 0], [1, 1] -1];       % RS 
%     [ [0, 0], [1, 0] -1]; % BS
%     [ [0, 1], [1, 1] -1];       % TS
];

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

c = 1; 

for j = 1:N
	for i = 1:N
        x0 = [ i/N j/N ];
        % solves a Laplace equation Δu = f at x0, where the boundary is given
        % by a collection of segments, and the boundary conditions are given
        % by a function g that can be evaluated at any point in space
        % screenedpoissonsolve(x0, segments, f, g, c)
        u = screenedpoissonsolve(x0, scene, @laplace_urefv, @newG, c);
        out(i,j) = u;
	end
end

disp("fin");
imagesc(out);

% New g function for input
% this one instead calculates the distance to a 'point' (any line
% segment), and returns 1 if we are close enough to it (epsilon), otherwise returns 0
function c = newG(xv, segments)
	%c = cos(2 * pi * x(1)) * sin(2 * pi * x(2)); 
    eps = 0.2;
    [h, w] = size(segments);
    
    distances = zeros(size(xv,1), h); 
    for j = 1:h
    	pv = closestpoints(xv, segments(j:j,:));
        distances(:,j:j) = vecnorm((xv-pv).').';
    end
    
    
    [M, ~] = min(distances,[],2);
    c = zeros(size(M,1),1);
    c(M(:) < eps) = 1;
    % lets just use 1 for now
    %s = segments(I(j),w);
    %c(j,1) = s; 
end


% reference solution (g)
function c = uref(xv, segments)
	%c = cos(2 * pi * x(1)) * sin(2 * pi * x(2)); 
   
    [h, w] = size(segments);
    
    distances = zeros(size(xv,1), h); 
    for j = 1:h
    	pv = closestpoints(xv, segments(j:j,:));
        distances(:,j:j) = vecnorm((xv-pv).').';
    end
    
    [~, I] = min(distances,[],2);
    
    c = zeros(size(I,1),1);
    for j = 1:size(I,1)
        % how do we want to treat the value of the boundary 
        s = segments(I(j),w);
        % just set it = to boundary condition?
        c(j,1) = s; 
        % trying out things
        % c(j,1) = fmod(floor(s .* xv(1)) + floor(s .* xv(2)), 2);
        % maybe we scale it by the distance from that boundary??? - no
        % c(j,1) = s * distances(j, I(j));
    end
end

% Laplacian of reference solution (f)
function c = laplace_urefv(xv)
    %c = 1; %8 .* (pi .* pi) .* cos(2.*pi.*xv(:,1:1)) .* sin(2 .* pi .* xv(:,2:2));
    c(xv(:,2) > 40) = 0;
    c(xv(:,2) <= 40) = 0;
    c = c';
end