% do some occupancy grid stuff
gridSize = 100; % 129 is the same as the WoSLaplace2D
n = gridSize*gridSize;


% Make a grid for occupancy of the profile
X = zeros(gridSize,gridSize);

% here is our scene (the walls which will be +ve)
lw = [[0.3, 0.3]; [0.3, 0.7]]*gridSize; % left wall
rw = [[0.7, 0.3]; [0.7, 0.7]]*gridSize; % right wall
tw = [[0.3, 0.3]; [0.7, 0.3]]*gridSize; % top wall
bw = [[0.45, 0.7]; [0.7, 0.7]]*gridSize; % bottom wall (has opening)
val = 1;    % set to 1

X = drawline(X, lw(1:1,:), lw(2:2,:), val);
X = drawline(X, rw(1:1,:), rw(2:2,:), val);
X = drawline(X, bw(1:1,:), bw(2:2,:), val);
X = drawline(X, tw(1:1,:), tw(2:2,:), val);

%imagesc(X);

% % lets draw lines around the outside that are -ve
% l = [[1, 1];[1,gridSize]];
% r = [[gridSize,1];[gridSize,gridSize]];
% b = [[1,1];[gridSize, 1]];
% t = [[1,gridSize];[gridSize,gridSize]];
% val = -1;   % set to -1
% 
% X = drawline(X, l(1:1,:), l(2:2,:), val);
% X = drawline(X, r(1:1,:), r(2:2,:), val);
% X = drawline(X, b(1:1,:), b(2:2,:), val);
% X = drawline(X, t(1:1,:), t(2:2,:), val);
% 
% imagesc(X);

rhsx = [];
C = [];
IX = @(i,j) (i-1)*gridSize+(j-1)+1;

for i = 1:gridSize    % row
    for j = 1:gridSize-1  % col
        if ( X(i,j) == 1 )
            rhsx = [rhsx; 1]; 
            crow = zeros(1,n);
            crow(1,IX(i,j)) = 1;
            C = [C; crow];      
        end
    end
end

subplot(1,2,1)
imagesc(X);
axis equal;
axis([1,gridSize,1,gridSize])

bdyConstraints = size(C,1); 
d = n + bdyConstraints;
S = sparse(d,d);

% http://12000.org/my_notes/mma_matlab_control/KERNEL/KEse82.htm
internalPoints=gridSize;
e   = ones(internalPoints,1);
spe = spdiags([e -2*e e], -1:1,internalPoints,internalPoints);
Iz  = speye(internalPoints);
S(1:n,1:n) = kron(Iz,spe)+kron(spe,Iz);

S(n+1:end,1:n) = C;
S(1:n,n+1:end) = C';

rhsxx = [ zeros(n,1); rhsx ];

% phix = conjgrad( S, rhsxx, 1e-8 );
phix = S\rhsxx;

norm(S*phix-rhsxx)

% reshape to draw in a grid
vx = reshape(phix(1:n),gridSize,gridSize)';

subplot(1,2,1)
%imagesc(vx);
%imagesc(log(1.00001-vx));
axis equal;
axis([1,gridSize,1,gridSize])

subplot(1,2,2)
imagesc(vx);
imagesc(log(1.000000001-vx).^2);
axis equal;
axis([1,gridSize,1,gridSize])