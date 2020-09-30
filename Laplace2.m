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

imagesc(X);

% lets draw lines around the outside that are -ve
l = [[0.01,0.01];[0.01,1]]*gridSize;
r = [[1,0.01];[1,1]]*gridSize;
b = [[0.01,0.01];[1,0.01]]*gridSize;
t = [[0.01,1];[1,1]]*gridSize;
val = -1;   % set to -1

X = drawline(X, l(1:1,:), l(2:2,:), val);
X = drawline(X, r(1:1,:), r(2:2,:), val);
X = drawline(X, b(1:1,:), b(2:2,:), val);
X = drawline(X, t(1:1,:), t(2:2,:), val);

imagesc(X);

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
