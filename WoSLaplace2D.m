%%%%%%%%%%%%%%%% data %%%%%%%%%%%%%%%%

% a segment is a pair of points
% the scene is an array(vector) of segments
scene = [
    [ [0.5, 0.1], [0.9, 0.5] ] ;
    [ [0.5, 0.9], [0.1, 0.5] ] ;
    [ [0.1, 0.5], [0.5, 0.1] ] ;
    [ [0.5, 0.33333333], [0.5, 0.6666666] ] ;
    [ [0.33333333, 0.5], [0.6666666, 0.5] ] ;
];
   
%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%
N = 129;    % image size
out = zeros(N,N);
out2 = zeros(N,N);

for j = 1:N
    fprintf("row %i of %i\n", j, N);
    for i = 1:N
        x0 = [ (i-1)./(N-1), (j-1)./(N-1) ];
        u = laplacesolve(x0, scene, @checker);
        out(i,j) = u;
    end
end

imagesc(out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = checker(x)
% boundary condition checker from the sample code 
    s = 6;
    c = fmod(floor(s .* x(1)) + floor(s .* x(2)), 2 );
end