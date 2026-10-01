function [F, info] = functionSweptMoment(P_ris, phat, Sigma, delta, vhat, lambda, r, degP, J)
%Low-rank factor of the second moment under BOTH prediction error and the
%displacement swept while one configuration is held (Theorem 2), built WITHOUT
%ever forming the N x N kernel.
%
%   Rbar = (a a^H) .* W,   W = W_pred .* W_sweep
%   W_pred_nm  = exp(-0.5 k0^2 (u_n-u_m)' Sigma (u_n-u_m))
%   W_sweep_nm = sinc( 0.5 k0 delta (u_n-u_m)' vhat ),   sinc(x) = sin(x)/x
%
%CONSTRUCTION. Each factor gets an explicit O(N r) feature map, and the Hadamard
%product of two low-rank matrices is factored by the columnwise (Khatri-Rao)
%product of their factors, since (F G^H) .* (P Q^H) = (F *kr* P)(G *kr* Q)^H.
%
%  Gaussian factor: the multi-index map of Theorem 3, dimension nchoosek(degP+3,3),
%  valid while t = 2 pi^2 zeta^2 stays small.
%
%  Sweep factor: sinc has the integral representation
%       sinc(c(g_n - g_m)) = (1/2) int_{-1}^{1} exp( j c (g_n - g_m) u ) du,
%  so a J-point Gauss-Legendre rule gives the exact-to-quadrature rank-J factor
%       F_s(n,j) = sqrt(w_j/2) * exp( j c g_n u_j ),      g_n = u_n' vhat,
%  converging exponentially in J for the arguments arising here. This is why the
%  sweep costs a small constant and not a dense eigendecomposition: an earlier
%  version formed W and ran a randomised range finder, which measured 585
%  iteration-equivalents of fixed cost and made the design look infeasible
%  everywhere on grounds that were an artefact of the implementation.
%
%The Khatri-Rao factor has rank r_p*J, most of which is redundant (the measured
%numerical rank of the product is at most two), so it is compressed to r by a thin
%QR -- O(N (r_p J)^2), still with no N x N matrix anywhere.
%
%This is version 2.0 (Last edited: 2026-09-11)
%License: GPLv2.

if nargin < 7 || isempty(r),    r    = 4; end
if nargin < 8 || isempty(degP), degP = 1; end
if nargin < 9 || isempty(J),    J    = 6; end

k0 = 2*pi/lambda;
a  = functionArrayResponse(P_ris, phat, lambda);
N  = numel(a);
V  = phat - P_ris;
U  = (V ./ vecnorm(V,2,1)).';
Ut = U - mean(U,1);

% ---- Gaussian factor, explicit feature map (Theorem 3)
[Qe,De] = eig((Sigma+Sigma.')/2,'vector');
Sh = Qe*diag(sqrt(max(real(De),0)))*Qe.';
Z  = k0*(Ut*Sh);
d  = exp(-0.5*sum(Z.^2,2));
Phi = ones(N,1);
for q = 1:degP
    for i = 0:q
        for jj = 0:(q-i)
            kk = q-i-jj;
            Phi = [Phi, (Z(:,1).^i).*(Z(:,2).^jj).*(Z(:,3).^kk) ...
                        /sqrt(factorial(i)*factorial(jj)*factorial(kk))];  %#ok<AGROW>
        end
    end
end
Fp = d .* Phi;                                   %N x r_p

% ---- sweep factor, Gauss-Legendre on the integral representation
if delta > 0
    g = Ut*vhat(:);                              %centring is harmless: only differences enter
    c = 0.5*k0*delta;
    [xg, wg] = localGauss(J);
    Fs = exp(1j*c*(g*xg.')) .* repmat(sqrt(wg.'/2), N, 1);   %N x J
else
    Fs = ones(N,1);
end

% ---- Khatri-Rao (columnwise) product, then compress
rp = size(Fp,2); rs = size(Fs,2);
Fkr = zeros(N, rp*rs);
c0 = 0;
for i = 1:rp
    Fkr(:, c0+(1:rs)) = Fp(:,i) .* Fs;
    c0 = c0 + rs;
end
Fkr = a .* Fkr;                                  %fold in the response

[Qk, Rk] = qr(Fkr, 0);
[Uk, Sk] = svd(Rk, 'econ');
rr = min(r, size(Sk,1));
F  = Qk*Uk(:,1:rr)*Sk(1:rr,1:rr);

if nargout > 1
    sv = diag(Sk);
    info = struct('rank', rr, 'rankFull', rp*rs, ...
                  'energyKept', sum(sv(1:rr).^2)/max(sum(sv.^2),eps), ...
                  't', max(sum(Z.^2,2)));
end
end

function [x,w] = localGauss(n)
%Gauss-Legendre nodes and weights on [-1,1] via the Golub-Welsch eigenproblem.
k = 1:n-1;
b = k./sqrt(4*k.^2-1);
[Vv,Dd] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(Dd));
w = 2*(Vv(1,ix).^2).';
x = x(:); w = w(:);
end
