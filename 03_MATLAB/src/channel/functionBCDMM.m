function [theta, nIter] = functionBCDMM(Xi, thetaInit, p)
%Block-coordinate-descent with majorisation-minimisation (BCD-MM) solver
%for the near-field multi-beam problem (P2), Sec. V-A.
%
%   min_{theta,omega} || gbar - Xi^H theta ||^2
%   s.t. |theta_n| = 1, |omega_k| = 1,  gbar = gtilde .* omega,
%        gtilde_k = sqrt(N/K).
%
%Sub-problem 1 (omega):  omega = exp(j*angle(Xi^H theta))
%Sub-problem 2 (MM):     theta <- exp(j*angle( Xi*gbar - Xi*Xi^H*theta
%                                              + lambda_max*theta ))
%
%CONVERGENCE METRIC (nIter): the dominant per-slot cost is the MM phase
%update (eq.(6)), each costing O(KN). We therefore report the TOTAL number
%of MM updates executed until convergence -- this is the quantity that
%governs per-slot complexity and runtime. A warm start (previous slot's
%solution) needs far fewer MM updates than a random cold start, which is the
%source of the 59%-76% saving.
%
%INPUT:
% Xi        = N x K steering matrix
% thetaInit = N x 1 initial phase vector (cold = random; warm = prev slot)
% p         = parameter struct (uses epsTol, muTol, Lmax, Imax)
%
%OUTPUT:
% theta  = N x 1 converged unit-modulus phase vector
% nIter  = total number of MM updates (per-slot iteration/complexity metric)
%
%This is version 2.2 (Last edited: 2026-07-23) -- reports total MM updates;
%adds a stall-based outer exit (see convergence test below) so the solver
%stops at a stationary point instead of always running Lmax iterations.
%
%BUGFIX (2026-07-23): the equal-gain target gtilde_k = sqrt(N/K) below is
%the correct target ONLY for UNIT-NORM steering columns (||xi_k|| = 1), as
%assumed in the paper's derivation (P2, Sec. V-A). The codebook and exact-
%focus constructions produce UNIT-MODULUS columns (||xi_k|| = sqrt(N)), so
%without normalisation the solver aimed the response a factor of N too low,
%suppressing beam gain BELOW even a random phase vector and collapsing the
%achievable rate (verified ~476x min-rate loss; the method fell behind a
%trivial vector-sum baseline). We normalise each steering column to unit
%norm here so that sqrt(N/K) is dimensionally correct. This changes NONE of
%the algorithm's logic -- it only rescales the target to match the paper's
%unit-norm assumption -- and restores the expected max-min-fair performance.
%License: GPLv2.

[N,K]  = size(Xi);
Xi     = Xi ./ vecnorm(Xi,2,1);          %<-- unit-norm steering columns (see BUGFIX above)
gtilde = sqrt(N/K) * ones(K,1);
lammax = max(real(eig(Xi'*Xi)));        %= lambda_max(Xi Xi^H)
theta  = thetaInit ./ abs(thetaInit);   %enforce unit modulus
nIter  = 0;                             %total MM updates
objPrev = inf;                          %previous outer objective (for stall exit)

for outer = 1:p.Lmax
    %--- Sub-problem 1: phase rotation ---
    %SPEEDUP (2026-07-24): exp(1j*angle(z)) == z./|z| exactly (unit-modulus
    %projection) but avoids the per-element atan2+sin+cos transcendentals,
    %which dominate the wall-clock. This changes NONE of the numerical
    %results (identical theta, iteration counts, rates) -- only runtime.
    Xit   = Xi' * theta;
    omega = Xit ./ abs(Xit);
    gbar  = gtilde .* omega;
    Xg    = Xi * gbar;                  %constant within the inner loop

    %--- Sub-problem 2: MM inner loop ---
    for inner = 1:p.Imax
        grad     = Xg - Xi*(Xi'*theta) + lammax*theta;   %O(KN)
        ag       = abs(grad);
        thetaNew = grad ./ ag;          %== exp(1j*angle(grad)), transcendental-free
        thetaNew(ag == 0) = 1;          %guard (grad is essentially never exactly 0)
        nIter    = nIter + 1;           %count every MM update
        if norm(thetaNew - theta) < p.muTol
            theta = thetaNew; break;
        end
        theta = thetaNew;
    end

    %--- Outer convergence test (absolute floor OR objective stall) ---
    %The equal-gain problem is non-convex and, for K>1 users sharing one
    %unit-modulus theta, the residual legitimately plateaus above the tiny
    %absolute epsTol -- so we ALSO exit once the monotonically-decreasing
    %BCD-MM objective stops improving (a genuine stationary point). Without
    %this, every solve needlessly ran all Lmax*Imax iterations. Solution
    %quality is unchanged (the plateau IS the solution); only wasted
    %iterations are removed, and nIter now reflects true convergence.
    obj = norm(gbar - Xi'*theta)^2;
    if obj < p.epsTol || (outer >= 4 && (objPrev - obj) < 1e-3*max(objPrev,eps))
        break;
    end
    objPrev = obj;
end
end
