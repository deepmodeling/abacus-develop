% GEN_BESSEL_GRID  Generate a Bessel function of the first kind in grid.
% q: the argument of the Bessel function
% rc: the cutoff radius
% nr: the number of points to be generated
% l_lo, l_hi: order of the Bessel function
%% Case: q = 0.1; rc = 50; nr = 5000; l_lo = 0; l_hi=11;
function [X, Val] = gen_bessel_grid(q, rc, nr, l_lo, l_hi)
    X = linspace(0, rc, nr+1) * q;
    X = X(2:nr+1);
    Val = zeros(nr, l_hi-l_lo+1);
    for l = l_lo:l_hi
        Val(:, l-l_lo+1) = sphbes(l, X);
    end
    Val = reshape(Val, 1, []);
end


