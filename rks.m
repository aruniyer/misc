function [Krr_approx, Z] = rks(X, Y, bw, n)

d = size(X, 2);
Z = randn(n, d)*(1 / bw);
phi_X = exp(i*Z*X') / sqrt(n);
phi_Y = exp(i*Z*Y') / sqrt(n);
Krr_approx = real(phi_X' * phi_Y);

end
