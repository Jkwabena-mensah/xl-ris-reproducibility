addPaper2Paths();
%exp12a_rate_outage_safe
%
%RESTART-SAFE RE-RUN OF exp04a (Section V-E, Figs 6 and 7).
%
%WHY. exp04a selects BOTH designs as the best of nRest = 2 restarts, and scores
%both on Rr -- the zeta-DEPENDENT moment. The nominal design does not depend on
%zeta, but its selection did, so the baseline configuration is free to move
%along the sweep. That is structurally the confound that manufactured the
%scattering reversal in Section V-G (see exp11b_v5g.m). This script removes it.
%
%THE FIX, and it is small. The nominal moment Rn does not depend on zeta, so the
%nominal designs are solved ONCE from a fixed restart bank and reused unchanged
%at every operating point. Nothing is selected. The comparison is then formed as
%a difference WITHIN each restart and WITHIN each evaluation draw, and reported
%as a mean over restarts with a standard error over restarts. There is no
%best-of anywhere in the protocol.
%
%WHAT TO EXPECT, stated before running so it cannot be rationalised afterwards.
%In exp04a the nominal design is handed its best restart AS JUDGED ON THE
%UNCERTAINTY-AWARE CRITERION -- the baseline is given every advantage. The bias
%therefore runs AGAINST the paper's claim, and the published gains should be
%understated rather than manufactured. If the restart-safe gains come out LOWER
%than published, the confound was working the other way and Section V-E needs
%the same treatment Section V-G got.
%
%PROTOCOL A is also computed, restricted to the first two restarts of the bank,
%to confirm that this script reproduces exp04a's published numbers before any
%weight is put on the Protocol B figures. Same discipline as exp11b.
%
%GEOMETRY. rng(11) and the same draw order as exp04a, so phat, rho and w_NF are
%identical; the restart bank is drawn immediately after phat and therefore does
%not disturb it.
%
%Version 1.0 (2026-09-18). License: GPLv2.

clear; close all; rng(11,'twister');
fprintf('\n===== exp12a: restart-safe re-run of exp04a =====\n');

p = functionSimParams(); P_ris = functionRISPositions(p);
N = p.N; K = p.K; lam = p.lambda;
pOpt = p; pOpt.Imax = 150; pOpt.muTol = 1e-4;

phat = zeros(3,K);
for k = 1:K
    phat(:,k) = [p.region(1,1)+diff(p.region(1,:))*rand;
                 p.region(2,1)+diff(p.region(2,:))*rand;
                 p.region(3,1)+diff(p.region(3,:))*rand];
end
rho = mean(vecnorm(phat,2,1));
wNF = lam*rho/p.D;

zetaG = [0.02 0.04 0.06 0.08 0.10 0.12];
nZ    = numel(zetaG);
nEval = 3000;
RB    = 12;                 % restart bank, matching exp11b
Rth   = 2;
SNRop = 10;
TH0   = exp(1j*2*pi*rand(N,RB));

fprintf('  N=%d K=%d rho=%.2f m w_NF=%.5f m | %d restarts, %d draws\n', ...
        N, K, rho, wNF, RB, nEval);

%% ---- nominal designs: solved ONCE, reused at every zeta ----------------
Rn = zeros(N,N,K);
for k = 1:K, Rn(:,:,k) = functionRobustMoment(P_ris, phat(:,k), zeros(3), lam); end
THn = zeros(N,RB);
tA = tic;
for r = 1:RB, THn(:,r) = localMaxMin(Rn, TH0(:,r), pOpt, N, K); end
fprintf('  nominal bank solved once (%.0f s)\n', toc(tA));

%% ---- sweep -------------------------------------------------------------
dRate = zeros(nZ,RB);        % paired min-rate difference, per restart
rNom  = zeros(nZ,RB); rAwa = zeros(nZ,RB);
oNom  = zeros(nZ,RB); oAwa = zeros(nZ,RB);
pA    = zeros(nZ,4);         % protocol A: [rNom rAwa oNom oAwa], best of 2
scale = NaN;

for iz = 1:nZ
    sig = zetaG(iz)*wNF;  Sig = sig^2*eye(3);
    Rr = zeros(N,N,K);
    for k = 1:K, Rr(:,:,k) = functionRobustMoment(P_ris, phat(:,k), Sig, lam); end

    THa = zeros(N,RB);
    for r = 1:RB, THa(:,r) = localMaxMin(Rr, TH0(:,r), pOpt, N, K); end

    TH = [THn THa];                                   % N x 2RB
    es = RandStream('twister','Seed',7000+iz);        % identical draws, all designs
    G  = zeros(nEval,K,2*RB);
    for m = 1:nEval
        At = functionArrayResponse(P_ris, phat + sig*randn(es,3,K), lam);
        G(m,:,:) = abs(At'*TH).^2/N^2;
    end

    %link budget fixed once, on the nominal bank at the smallest zeta
    if isnan(scale)
        gRef  = median(reshape(G(:,:,1:RB),[],1));
        scale = 10^(SNRop/10)/gRef;
        fprintf('  gRef = %.6e (restart-averaged nominal at zeta = %.2f)\n', gRef, zetaG(1));
    end

    Rm = log2(1 + scale*G);                           % nEval x K x 2RB
    mn = squeeze(min(Rm,[],2));                       % nEval x 2RB
    rNom(iz,:) = mean(mn(:,1:RB),1);
    rAwa(iz,:) = mean(mn(:,RB+1:end),1);
    %paired at the draw level as well as the restart level
    dRate(iz,:) = mean(mn(:,RB+1:end) - mn(:,1:RB), 1);
    ob = squeeze(mean(mean(Rm < Rth,1),2))';          % 1 x 2RB
    oNom(iz,:) = ob(1:RB); oAwa(iz,:) = ob(RB+1:end);

    %protocol A: best of the first two restarts, scored on Rr -- exp04a's rule
    bn = -inf; ba = -inf; ia = 1; ib = 1;
    for r = 1:2
        en = inf; ea = inf;
        for k = 1:K
            en = min(en, real(THn(:,r)'*Rr(:,:,k)*THn(:,r))/N^2);
            ea = min(ea, real(THa(:,r)'*Rr(:,:,k)*THa(:,r))/N^2);
        end
        if en > bn, bn = en; ia = r; end
        if ea > ba, ba = ea; ib = r; end
    end
    pA(iz,:) = [rNom(iz,ia) rAwa(iz,ib) oNom(iz,ia) oAwa(iz,ib)];

    fprintf('  zeta %.2f | B: nom %.4f awa %.4f d %+.4f | A: nom %.4f awa %.4f | %.0f s\n', ...
        zetaG(iz), mean(rNom(iz,:)), mean(rAwa(iz,:)), mean(dRate(iz,:)), ...
        pA(iz,1), pA(iz,2), toc(tA));
    save('exp12a_safe.mat','zetaG','dRate','rNom','rAwa','oNom','oAwa','pA', ...
         'scale','RB','nEval','Rth','rho','wNF','phat','iz');
end

%% ---- report ------------------------------------------------------------
se = @(x) 1.96*std(x,0,2)/sqrt(RB);
sd = se(dRate); md = mean(dRate,2); td = md./(sd/1.96);
dO = oAwa - oNom; mo = mean(dO,2); so = se(dO);

fprintf('\n PROTOCOL B (restart-safe): paired min-rate gain, mean +- 95%% CI over %d restarts\n', RB);
fprintf(' zeta |  nominal   aware   |    gain   +- CI      t   | outage nom  awa   paired d\n');
for iz = 1:nZ
    fprintf(' %.2f | %7.4f %7.4f | %+7.4f %7.4f %6.2f | %9.3e %9.3e %+9.3e\n', ...
        zetaG(iz), mean(rNom(iz,:)), mean(rAwa(iz,:)), md(iz), sd(iz), td(iz), ...
        mean(oNom(iz,:)), mean(oAwa(iz,:)), mo(iz));
end
fprintf('\n resolved above zero at %d of %d operating points (CI excludes 0)\n', ...
        sum(md-sd > 0), nZ);
fprintf(' outage reduction resolved at %d of %d\n', sum(mo+so < 0), nZ);

fprintf('\n PROTOCOL A (exp04a''s rule, best of 2 scored on Rr) -- the reproduction check\n');
fprintf(' zeta |  nominal   aware  |    gain\n');
for iz = 1:nZ
    fprintf(' %.2f | %7.4f %7.4f | %+7.4f\n', zetaG(iz), pA(iz,1), pA(iz,2), pA(iz,2)-pA(iz,1));
end

pub = load('exp04a_rate_outage.mat');
fprintf('\n PUBLISHED (exp04a_rate_outage.mat)\n');
fprintf(' zeta |  nominal   aware  |    gain\n');
for iz = 1:nZ
    fprintf(' %.2f | %7.4f %7.4f | %+7.4f\n', ...
        pub.zetaG(iz), pub.minR(iz,1), pub.minR(iz,2), pub.minR(iz,2)-pub.minR(iz,1));
end

done = true;
save('exp12a_safe.mat','zetaG','dRate','rNom','rAwa','oNom','oAwa','pA', ...
     'scale','RB','nEval','Rth','rho','wNF','phat','done');
fprintf('\n===== saved exp12a_safe.mat =====\n');

%% ================= local functions =================
function th = localMaxMin(R, th0, pOpt, N, K)
%Byte-identical in behaviour to exp04a's local solver: max-min by bisection on
%the common target, the inner equal-gain iteration used only as a feasibility
%test at a fixed target (F-35).
  lo = sqrt(N/K); hi = N; th = functionRobustEqualGainMM(R, th0, lo, pOpt);
  for it = 1:9
    gb = 0.5*(lo+hi);
    t  = functionRobustEqualGainMM(R, th0, gb, pOpt);
    psi = zeros(K,1);
    for k = 1:K, psi(k) = sqrt(max(real(t'*R(:,:,k)*t),0)); end
    if min(psi) >= 0.98*gb, lo = gb; th = t; else, hi = gb; end
  end
end
