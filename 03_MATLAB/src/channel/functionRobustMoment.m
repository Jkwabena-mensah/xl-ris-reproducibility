function [R, Rfac, W] = functionRobustMoment(P_ris, phat, Sigma, lambda, rmax)
%Closed-form second moment of the near-field array response under Gaussian
%position uncertainty:  R = E[ a(phat+e) a(phat+e)^H ],  e ~ N(0, Sigma).
%
%Using the linearised response of Lemma 1, a(phat+e) = a(phat) .* exp(-j k0 U e),
%the Gaussian characteristic function gives, exactly at that order,
%
%   R_nm = a_n a_m^*  exp( -0.5 k0^2 (u_n - u_m)' Sigma (u_n - u_m) )
%
%so R = (a a^H) .* W with W real, symmetric and positive definite. No sampling
%is required: the sample count S disappears from the design, and with it the
%Monte-Carlo error in the moment itself.
%
%Sigma -> 0 gives W -> ones(N) and R -> a a^H, so every quantity built on R
%reduces to its nominal counterpart. This reduction is the primary validation.
%
%INPUT:
% P_ris  = 3 x N RIS element positions [m]
% phat   = 3 x 1 predicted user position [m]
% Sigma  = 3 x 3 position-error covariance [m^2]
% lambda = wavelength [m]
% rmax   = (optional) rank for the truncated factor Rfac; [] or omitted
%          returns Rfac = [] and skips the eigendecomposition
%
%OUTPUT:
% R    = N x N Hermitian PSD second-moment matrix
% Rfac = N x rmax factor with Rfac*Rfac' ~= R (empty if rmax not supplied).
%        R is numerically rank <= 3 for zeta <= 0.12, so a small rmax loses
%        nothing; this is what keeps the per-iteration cost O(rN) rather
%        than O(N^2).
% W    = N x N real decoherence kernel
%
%This is version 1.0 (Last edited: 2026-09-07)
%License: GPLv2.

if nargin < 5, rmax = []; end

k0 = 2*pi/lambda;
a  = functionArrayResponse(P_ris, phat, lambda);      %N x 1

V = phat - P_ris;                                     %3 x N
U = (V ./ vecnorm(V,2,1)).';                          %N x 3, rows u_n'

%(u_n - u_m)' Sigma (u_n - u_m) = q_n + q_m - 2 u_n' Sigma u_m
q     = sum((U*Sigma).*U, 2);                         %N x 1
Cross = U*Sigma*U.';
Dq    = q + q.' - 2*Cross;
Dq    = max(Dq, 0);                                   %guard tiny negatives

W = exp(-0.5*k0^2*Dq);
W = (W + W.')/2;

R = (a*a') .* W;
R = (R + R')/2;

Rfac = [];
if ~isempty(rmax)
    [Vv, Dd] = eig(R, 'vector');
    [Dd, ix] = sort(real(Dd), 'descend');
    Vv = Vv(:, ix);
    r  = min(rmax, sum(Dd > 0));
    Rfac = Vv(:,1:r) * diag(sqrt(max(Dd(1:r),0)));
end
end
