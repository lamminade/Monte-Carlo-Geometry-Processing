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
    [ [0.2, 0.2], [0.2, 0.7] 1]; % left wall (LW)
    [ [0.8, 0.2], [0.8, 0.7] 1]; % right wall
    [ [0.2, 0.2], [0.8, 0.2] 1]; % bottom wall
    [ [0.4, 0.7], [0.8, 0.7] 1]; % top wall (has opening)
     % outside perimeter of image ! 
    [ [0, 0], [0, 1] -1];    % LS
    [ [1, 0], [1, 1] -1];       % RS 
    [ [0, 0], [1, 0] -1]; % BS
    [ [0, 1], [1, 1] -1];       % TS
];

%%%%%%%%%%%%%%%% main - vector version %%%%%%%%%%%%%%%%
N = 129;           % image size
out = zeros(N,N);

for j = 1:N
    fprintf("row %i of %i\n", j, N);
    for i = 1:N
        x0 = [ (i-1)./(N-1), (j-1)./(N-1) ];
        u = laplacesolve(x0, scene, @checkerTwo);
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

% pls vectorize this
% function c = checkerTwoSlow(segments, xv)
% % new boundary condition checker where boundary condition is instead
% %   attached to the segments as the final value of the segment
% % checks which segment/edge is closest to x and then returns value of the
% %   bounding function specified on that edge
%     distances = zeros(1, size(segments, 1)); 
%     for i = 1:size(xv,1)
%         for j = 1:size(segments,1)
%             p = closestpoint(xv(i:i,:), segments(j:j,:));
%             distances(j) = vecnorm((xv(i:i,:)-p(i:i,:)).').'; %pdist([x;p],'euclidean') ?;
%         end
%     end
%     % what to do when 2 equally close segments ??
%     % for now .. i will just take the first -- bring this up w paul?
%     [~, i] = min(distances);
%     c = segments(i, size(segments,2)) * distances(i);
% end

% vectorized it
function c = checkerTwo(segments, xv)
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
        % whatever this is
        % c(j,1) = fmod(floor(s .* xv(1)) + floor(s .* xv(2)), 2);
        % maybe we scale it by the distance from that boundary???
        % c(j,1) = s * distances(j, I(j));
    end
end