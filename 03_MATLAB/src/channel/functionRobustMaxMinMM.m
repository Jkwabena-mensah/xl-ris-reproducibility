function [theta, nIter, hist] = functionRobustMaxMinMM(Rk, thetaInit, p)
%Uncertainty-aware max-min beam design for a passive XL-RIS:
%
%   max_theta  min_k  theta^H Rbar_k theta      s.t.  |theta_n| = 1,
%
%where Rbar_k = E[a(p_k) a(p_k)^H] is the closed-form second moment returned by
%functionRobustMoment. The objective is the EXPECTED focusing gain of user k,
%so it is invariant to any phase common to all elements -- which is precisely
%the defect of the sample-averaged minimum-distance surrogate, whose linear
%term -2Re{theta^H Xibar g} is annihilated by common-mode decoherence at any
%realistic sigma_p (the common-mode phase variance k0^2 ubar'Sigma ubar exceeds
%the aperture-spread variance k0^2 tr(C Sigma) by three orders of magnitude).
%
%ALGORITHM. Each f_k(theta) = theta^H Rbar_k theta is convex, hence minorised
%at theta^(t) by its linearisation
%
%   f_k(theta) >= 2 Re{theta^H Rbar_k theta^(t)} - theta^(t)H Rbar_k theta^(t).
%
%The non-smooth min over k is handled by the smooth surrogate
%F_mu = -mu log sum_k exp(-f_k/mu), whose gradient weights are the softmin
%weights w_k. Maximising the weighted minorant over the unit-modulus set has
%the closed form
%
%   theta^(t+1) = exp( j angle( sum_k w_k Rbar_k theta^(t) ) ),
%
%a weighted power iteration. mu is annealed downward so the surrogate tightens
%onto the true max-min. With Sigma -> 0 every Rbar_k -> a_k a_k^H and the
%update reduces to the nominal max-min focusing iteration.
%
%INPUT:
% Rk        = N x N x K array of second moments, OR 1 x K cell of N x r factors
%             Rfac_k (then Rbar_k*theta is evaluated as Rfac*(Rfac'*theta))
% thetaInit = N x 1 unit-modulus initialisation (warm start or random)
% p         = parameter struct (uses Imax, muTol)
%
%OUTPUT:
% theta  = N x 1 unit-modulus design
% nIter  = number of MM updates executed (the per-slot complexity metric,
%          consistent with functionBCDMM)
% hist   = struct with .minGain per iteration, for convergence plots
%
%This is version 1.0 (Last edited: 2026-09-07)
%License: GPLv2.

useFac = iscell(Rk);
if useFac, K = numel(Rk); N = size(Rk{1},1);
else,      K = size(Rk,3); N = size(Rk,1); end

theta = thetaInit(:);
theta = theta./abs(theta);

Imax = p.Imax;
if isfield(p,'muTol'), tol = p.muTol; else, tol = 1e-3; end

%Annealing schedule for the softmin temperature, relative to the gain scale.
mu0 = 0.10;  muEnd = 1e-3;

minGain = zeros(Imax,1);
nIter   = 0;

Rtheta = zeros(N,K);
for t = 1:Imax
    %---- per-user gains and Rbar_k * theta, reused by both steps
    f = zeros(K,1);
    for k = 1:K
        if useFac, Rtheta(:,k) = Rk{k}*(Rk{k}'*theta);
        else,      Rtheta(:,k) = Rk(:,:,k)*theta; end
        f(k) = real(theta'*Rtheta(:,k));
    end
    minGain(t) = min(f)/N^2;

    %---- softmin weights (shift-stabilised)
    mu = mu0*(muEnd/mu0)^((t-1)/max(Imax-1,1));
    sc = max(mean(f), eps);
    z  = -(f - min(f))/(mu*sc);
    w  = exp(z - max(z));
    w  = w/sum(w);

    %---- weighted minorant maximiser over the unit-modulus set
    grad     = Rtheta*w;
    thetaNew = exp(1j*angle(grad));
    nIter    = nIter + 1;

    if norm(thetaNew - theta) < tol && t > 5
        theta = thetaNew;
        minGain(t+1:end) = minGain(t);
        break
    end
    theta = thetaNew;
end

hist.minGain = minGain(1:max(nIter,1));
end
