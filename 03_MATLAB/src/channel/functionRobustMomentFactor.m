function [F, info] = functionRobustMomentFactor(P_ris, phat, Sigma, lambda, deg)
%Low-rank factor of the closed-form second moment, built in O(N r) WITHOUT
%forming the N x N matrix and WITHOUT an eigendecomposition.
%
%   Rbar = E[a(phat+e) a(phat+e)^H] ~= F F^H,   F in C^{N x r}.
%
%DERIVATION. From Theorem 2, Rbar = (a a^H) .* W with
%W_nm = exp(-0.5 k0^2 (u_n-u_m)' Sigma (u_n-u_m)). The mean direction cancels in
%the difference, so with utilde_n = u_n - ubar and z_n = k0 Sigma^{1/2} utilde_n,
%
%   W_nm = d_n d_m exp(z_n . z_m),      d_n = exp(-0.5 ||z_n||^2).
%
%The z_n are confined to a cone of half-angle ~D/2rho, and
%
%   max_n ||z_n||^2 = 2 pi^2 zeta^2                                        (*)
%
%exactly, where zeta = sigma_p/wNF -- the SAME dimensionless group, and the same
%factor 2 pi^2, that sets the focusing loss in Corollary 1. Because that argument
%is small over the whole validity range (t = 0.0079 at zeta = 0.02, t = 0.284 at
%zeta = 0.12), the exponential truncates after very few terms:
%
%   exp(z.z') = sum_{q<=p} (z.z')^q/q! + E,   |E| <= sum_{q>p} t^q/q!,  t as in (*).
%
%Truncating at degree p gives an EXACT feature map of dimension nchoosek(p+3,3)
%-- 4 for p=1, 10 for p=2 -- so the rank bound is derived rather than measured.
%Measured relative Frobenius error at p=2: 6.1e-9 (zeta=0.02), 2.4e-5 (0.08),
%2.6e-4 (0.12), against the bound 8.2e-8, 3.5e-4, 4.1e-3.
%
%The feature for multi-index alpha is z^alpha / sqrt(alpha!), which reproduces
%(z.z')^q/q! term by term.
%
%INPUT:
% P_ris  = 3 x N RIS element positions [m]
% phat   = 3 x 1 predicted user position [m]
% Sigma  = 3 x 3 position-error covariance [m^2]
% lambda = wavelength [m]
% deg    = (optional) truncation degree p, default 2
%
%OUTPUT:
% F    = N x r factor, r = nchoosek(deg+3,3), with F*F' ~= Rbar
% info = struct: .t (= 2 pi^2 zeta^2, the series argument), .rank, .tailBound
%
%This is version 1.0 (Last edited: 2026-09-07)
%License: GPLv2.

if nargin < 5 || isempty(deg), deg = 2; end

k0 = 2*pi/lambda;
a  = functionArrayResponse(P_ris, phat, lambda);      %N x 1
N  = numel(a);

V  = phat - P_ris;                                    %3 x N
U  = (V ./ vecnorm(V,2,1)).';                         %N x 3
Ut = U - mean(U,1);

%Symmetric square root; Sigma is a covariance, hence PSD.
[Q,D] = eig((Sigma+Sigma.')/2, 'vector');
Sh    = Q*diag(sqrt(max(real(D),0)))*Q.';
Z     = k0*(Ut*Sh);                                   %N x 3
nz    = sum(Z.^2, 2);
d     = exp(-0.5*nz);

%Feature map: all monomials z^alpha/sqrt(alpha!) with |alpha| <= deg.
Phi = ones(N,1);
for q = 1:deg
    for i = 0:q
        for j = 0:(q-i)
            k = q - i - j;
            col = (Z(:,1).^i).*(Z(:,2).^j).*(Z(:,3).^k) ...
                  / sqrt(factorial(i)*factorial(j)*factorial(k));
            Phi = [Phi, col];                                        %#ok<AGROW>
        end
    end
end

F = (a .* d) .* Phi;                                  %N x r

%Guard: the truncation is only justified while the series argument is small.
%t = 2 pi^2 zeta^2, so t > 0.3 corresponds to zeta > 0.123 -- outside the stated
%validity range, where the reinitialisation criterion should already have fired.
tArg = max(nz);
tailB = exp(tArg) - sum(tArg.^(0:deg)./factorial(0:deg));
if tailB > 1e-2
    warning('functionRobustMomentFactor:truncation', ...
        ['degree %d leaves a relative error bound of %.2e (t = %.3f, zeta = %.3f). '...
         'Raise deg or reinitialise.'], deg, tailB, tArg, sqrt(tArg/(2*pi^2)));
end

if nargout > 1
    info = struct('t', tArg, 'rank', size(F,2), 'tailBound', tailB);
end
end
