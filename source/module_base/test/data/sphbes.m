function js = sphbes(nu, x)
% returns the spherical Bessel functions jnu(x)
% x is a vector or it may be a matrix if nu is a scalar
% if nu is a row and x a column vector, the output js is a matrix

[~, lnu] = size(nu);
xm = repmat(x, 1, lnu);
% special case handle: x = 0
js(~xm) = 0;
% general case
nzind = find(xm);
js(nzind) = sqrt(pi ./(2* xm(nzind))) .* besselj(nu + 0.5, x(nzind));




