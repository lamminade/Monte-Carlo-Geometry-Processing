function GrR = harmonicgreensv(r, R)
% HARMONICGREENSV - vectorized harmonic Green's function 
% for a 2D ball of radius R
	GrR = log(R./r) / (2 * pi);
    GrR(GrR == -inf) = 0;
end