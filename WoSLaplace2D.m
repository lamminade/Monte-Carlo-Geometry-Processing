%%%%%%%%%%%%%%%% data %%%%%%%%%%%%%%%%

% the original scene
% scene = [
%     [ [0.5, 0.1], [0.9, 0.5] ] ;
%     [ [0.5, 0.9], [0.1, 0.5] ] ;
%     [ [0.1, 0.5], [0.5, 0.1] ] ;
%     [ [0.5, 0.33333333], [0.5, 0.6666666] ] ;
%     [ [0.33333333, 0.5], [0.6666666, 0.5] ] ;
% ];

% Badpuzzle with boundary conditions attached
scene = [
    [ [0.3, 0.3], [0.3, 0.7] 1]; % top wall
    [ [0.7, 0.45], [0.7, 0.7] 1]; % bottom wall (has opening)
    [ [0.3, 0.3], [0.7, 0.3] 1]; % left wall
    [ [0.3, 0.7], [0.7, 0.7] 1]; % right wall 
     % outside perimeter of image ! 
    [ [0, 0], [0, 1] -1];    % LS
    [ [1, 0], [1, 1] -1];       % RS 
    [ [0, 0], [1, 0] -1]; % BS
    [ [0, 1], [1, 1] -1];       % TS
];

%%%%%%%%%%%%%%%% main - vector version %%%%%%%%%%%%%%%%
N = 100;           % image size
out = zeros(N,N);

for j = 1:N
    %fprintf("row %i of %i\n", j, N);
    for i = 1:N
        x0 = [ (i-1)./(N-1), (j-1)./(N-1) ];
        u = laplacesolve(x0, scene, @getboundaryvalue);
        out(i,j) = u;
    end
end

imagesc(out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boundary condition functions to pass into the solve

function c = checker(x)
% checker from the sample code 
    s = 6;
    c = fmod(floor(s .* x(1)) + floor(s .* x(2)), 2 );
end

function c = getboundaryvalue(segments, xv)
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