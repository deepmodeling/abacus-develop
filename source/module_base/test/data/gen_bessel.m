% GEN_BESSEL  Generate a Bessel function of the first kind.
% X: function domain of the Bessel function
% l_lo, l_hi: euqation order of the Bessel function
%% Case: X = 2.^[-5:-1:-20]; l_lo = 0; l_hi=11;
function Val = gen_bessel(X, l_lo, l_hi)
    [~, n] = size(X);
    Val = zeros(n, l_hi-l_lo+1);
    for l = l_lo:l_hi
        Val(:, l-l_lo+1) = sphbes(l, X);
    end
    Val = reshape(Val, 1, []);
end

