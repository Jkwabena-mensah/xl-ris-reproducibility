function [theta, nIter, hist] = functionRobustEqualGainMM(Rk, thetaInit, gbar, p)
%Uncertainty-aware equal-gain beam design for a passive XL-RIS.
%
%   min_theta  sum_k ( psi_k(theta) - gbar )^2 ,   psi_k = sqrt(theta^H Rbar_k theta),
%   s.t. |theta_n| = 1,
%
%where Rbar_k = E[a(p_k) a(p_k)^H] is the closed-form second moment from
%functionRobustMoment and psi_k is the root-mean-square expected amplitude
%delivered to user k. Driving every psi_k to the common target gbar = sqrt(N/K)
%is the fairness surrogate of Paper 1, carried over to the uncertain case.
%
%WHY NOT max-min DIRECTLY. A first attempt maximised min_k theta^H Rbar_k theta
%with softmin weights. It fails: as the softmin temperature is annealed the
%weights concentrate on the single worst user, the update focuses the whole
%surface on that user, the identity of the worst user then changes, and the
%iteration oscillates. Measured min-gain fell from 1.05e-3 to ~0 and the ascent
%was non-monotone (worst step -5.9e-3). The equal-gain surrogate below has no
%such failure mode because every user contributes at every iteration.
%
%WHY NOT THE SAMPLE-AVERAGED MINIMUM-DISTANCE FORM. Averaging
%||g - Xi^H theta||^2 over position draws with omega fixed outside the
%expectation is degenerate: its linear term carries Xibar = E[a], which is
%annihilated by common-mode phase decoherence (the common-mode phase variance
%k0^2 ubar' Sigma ubar exceeds the aperture-spread variance k0^2 tr(C Sigma) by
%three orders of magnitude), leaving an objective minimised by defocusing. The
%formulation here is invariant to any phase common to all elements, so that
%failure cannot occur.
%
%MM DERIVATION. Write J = theta^H R theta - 2 gbar sum_k psi_k + K gbar^2 with
%R = sum_k Rbar_k. The first term is convex and is majorised at theta^(t) by
%its lambdabar-quadratic upper bound; psi_k is convex (a norm of a linear map)
%and is therefore minorised by Re{theta^H Rbar_k theta^(t)}/psi_k^(t), which
%majorises -2 gbar psi_k. Minimising the resulting linear majorant over the
%unit-modulus set is closed form:
%
%   theta^(t+1) = exp( j angle( lambdabar theta^(t) - R theta^(t)
%                               + gbar sum_k Rbar_k theta^(t)/psi_k^(t) ) ).
%
%REDUCTION. At Sigma = 0, Rbar_k = a_k a_k^H and psi_k = |a_k^H theta|, so
%Rbar_k theta / psi_k = a_k exp(j angle(a_k^H theta)) and the update becomes
%exactly the BCD-MM iteration of Paper 1 with omega_k = exp(j angle(xi_k^H
%theta)). The auxiliary phase of Paper 1 is thus recovered rather than assumed.
%
%INPUT:
% Rk        = N x N x K array of second moments, OR 1 x K cell of N x r factors
% thetaInit = N x 1 unit-modulus initialisation (warm start or random)
% gbar      = common amplitude target, typically sqrt(N/K)
% p         = parameter struct (uses Imax, muTol)
%
%OUTPUT:
% theta = N x 1 unit-modulus design
% nIter = number of MM updates (per-slot complexity metric)
% hist  = struct with .obj and .minGain per iteration
%
%This is version 1.0 (Last edited: 2026-09-07)
%License: GPLv2.

useFac = iscell(Rk);
if useFac, K = numel(Rk); N = size(Rk{1},1);
else,      K = size(Rk,3); N = size(Rk,1); end

theta = thetaInit(:);  theta = theta./abs(theta);
Imax = p.Imax;
if isfield(p,'muTol'), tol = p.muTol; else, tol = 1e-4; end

%lambdabar >= lambda_max(sum_k Rbar_k). Gershgorin avoids an eigendecomposition
%and is valid as a majorisation constant (a loose bound only slows convergence).
if useFac
    Rsum = zeros(N);
    for k = 1:K, Rsum = Rsum + Rk{k}*Rk{k}'; end
else
    Rsum = sum(Rk,3);
end
lamBar = max(sum(abs(Rsum),2));

obj = zeros(Imax,1);  mg = zeros(Imax,1);  nIter = 0;
Rt = zeros(N,K);

for t = 1:Imax
    psi = zeros(K,1);
    for k = 1:K
        if useFac, Rt(:,k) = Rk{k}*(Rk{k}'*theta);
        else,      Rt(:,k) = Rk(:,:,k)*theta; end
        psi(k) = sqrt(max(real(theta'*Rt(:,k)), 0));
    end
    obj(t) = sum((psi - gbar).^2);
    mg(t)  = min(psi).^2/N^2;

    %majorant minimiser
    v = lamBar*theta - Rsum*theta + gbar*(Rt*(1./max(psi,eps)));
    thetaNew = exp(1j*angle(v));
    nIter = nIter + 1;

    if norm(thetaNew - theta) < tol && t > 5
        theta = thetaNew;
        obj(t+1:end) = obj(t);  mg(t+1:end) = mg(t);
        break
    end
    theta = thetaNew;
end

hist.obj     = obj(1:max(nIter,1));
hist.minGain = mg(1:max(nIter,1));
end
