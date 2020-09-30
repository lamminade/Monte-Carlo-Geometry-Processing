function [fin] = drawline(mat, p1, p2, val)
%DRAWLINE Draws a line from point p1 to point p2 on the matrix mat
%   Assumes the matrix mat is (mostly) zeros, draws a 'line' of pixels by 
% from point p1 to point p2 by setting them to val and returns a new matrix
% Adapted from https://stackoverflow.com/a/1941310
    fin = mat;  % idk how mutability works in matlab oop

    % calculate num pixels in line
    x = [p1(1) p2(1)];  % the x coordinates
    y = [p1(2) p2(2)];  % the y coordinates
    nPoints = max(abs(diff(x)), abs(diff(y)))+1;  % Number of points in line

    % Next, we compute row and column indices for the line pixels using 
    %   linspace, convert them from subscripted indices to linear indices 
    %   using sub2ind, then use them to modify mat:
    rIndex = round(linspace(y(1), y(2), nPoints));  % Row indices
    cIndex = round(linspace(x(1), x(2), nPoints));  % Column indices
    index = sub2ind(size(fin), rIndex, cIndex);     % Linear indices
    fin(index) = val;  % Set the line pixels to whatever u want (we want 1)
end

