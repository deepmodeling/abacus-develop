% SPHBES: Computes the Spherical Bessel functions of the first kind
% nu: euqation order, specified as a scalar, vector, matrix, or multidimensional array.
% x: function domain, specified as a scalar, vector, matrix, or multidimensional array.
function js = sphbes(nu, x)
    % special case handle: x = 0
    js(~x) = 0;
    % general case
    nzind = find(x);
    js(nzind) = sqrt(pi ./(2* xm(nzind))) .* besselj(nu + 0.5, x(nzind));
end




