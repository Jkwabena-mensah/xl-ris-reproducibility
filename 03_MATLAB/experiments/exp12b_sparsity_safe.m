addPaper2Paths();
%exp12b_sparsity_safe
%
%RESTART-SAFE RE-RUN OF exp04b (Section V-F, Fig. 8).
%
%WHY. exp04b selects both designs as the best of nRest = 5 restarts, scored on
%Rr -- the D_fb-DEPENDENT moment. The nominal design does not depend on D_fb but
%its selection did, which is the confound identified in exp11b_v5g.m. exp04b at
%least holds the restart BANK common across the sweep (TH0), which exp04a does
%not; the selection within that bank is still free to move.
%
%THE FIX. The nominal moment Rn does not depend on D_fb, so the nominal designs
%are solved ONCE from a fixed bank and reused unchanged at every operating
%point. The comparison is a difference within each restart and within each
%evaluation draw, averaged over restarts with a standard error over restarts.
%No best-of anywhere.
%
%DIRECTION OF THE EXPECTED BIAS, stated before running. exp04b hands the nominal
%design its best restart as judged on the uncertainty-aware criterion, so the
%baseline is favoured and the published gains should be understated. A LOWER
%restart-safe gain would mean the bias ran the other way and Section V-F needs
%correcting rather than confirming.
%
%GEOMETRY. rng(23) and the same draw order as exp04b, so phat, rho and w_NF are
%identical; the restart bank is drawn immediately after phat.
%
%Version 1.0 (2026-09-18). License: GPLv2.

clear; close all; rng(23,'twister');
fprintf('\n===== exp12b: restart-safe re-run of exp04b =====\n');

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

Ts   = 1e-3;
v    = 60/3.6;
vhat = [0;1;0];
DfbG = [1 2 5 10 20 50];
nD   = numel(DfbG);
nEval = 3000; RB = 12; Rth = 2; SNRop = 10;
TH0  = exp(1j*2*pi*rand(N,RB));

fprintf('  N=%d K=%d Ts=%.1f ms v=%.1f km/h w_NF=%.5f m | %d restarts, %d draws\n', ...
        N, K, Ts*1e3, v*3.6, wNF, RB, nEval);

%% ---- nominal designs: solved ONCE ---------------------------------------
Rn = zeros(N,N,K);
for k = 1:K, Rn(:,:,k) = functionRobustMoment(P_ris, phat(:,k), zeros(3), lam); end
THn = zeros(N,RB);
tA = tic;
for r = 1:RB, THn(:,r) = localMaxMin(Rn, TH0(:,r), pOpt, N, K); end
fprintf('  nominal bank solved once (%.0f s)\n', toc(tA));

%% ---- sweep --------------------------------------------------------------
dRate = zeros(nD,RB); rNom = zeros(nD,RB); rAwa = zeros(nD,RB);
oNom  = zeros(nD,RB); oAwa = zeros(nD,RB);
rmsePos = zeros(nD,1); zPred = zeros(nD,1); zSweep = zeros(nD,1);
pA = zeros(nD,4);
scale = NaN;

for id = 1:nD
    Dfb   = DfbG(id);
    Sig   = localCov(Ts, Dfb, 1.0, p.r);
    sigp  = sqrt(trace(Sig)/3);
    delta = v*Ts*Dfb;
    rmsePos(id) = sqrt(trace(Sig));
    zPred(id)   = sigp/wNF;
    zSweep(id)  = delta/wNF;

    Rr = zeros(N,N,K);
    for k = 1:K
        Fk = functionSweptMoment(P_ris, phat(:,k), Sig, delta, vhat, lam, 4, 2, 6);
        Rr(:,:,k) = Fk*Fk';
    end

    THa = zeros(N,RB);
    for r = 1:RB, THa(:,r) = localMaxMin(Rr, TH0(:,r), pOpt, N, K); end

    TH = [THn THa];
    es = RandStream('twister','Seed',4242);       % identical draws at every D_fb
    G  = zeros(nEval,K,2*RB);
    for m = 1:nEval
        e  = sigp*randn(es,3,K) + vhat*((rand(es,1,K)-0.5)*delta);
        At = functionArrayResponse(P_ris, phat + e, lam);
        G(m,:,:) = abs(At'*TH).^2/N^2;
    end

    if isnan(scale)
        gRef  = median(reshape(G(:,:,1:RB),[],1));
        scale = 10^(SNRop/10)/gRef;
        fprintf('  gRef = %.6e (restart-averaged nominal at D_fb = %d)\n', gRef, DfbG(1));
    end

    Rm = log2(1 + scale*G);
    mn = squeeze(min(Rm,[],2));
    rNom(id,:)  = mean(mn(:,1:RB),1);
    rAwa(id,:)  = mean(mn(:,RB+1:end),1);
    dRate(id,:) = mean(mn(:,RB+1:end) - mn(:,1:RB), 1);
    ob = squeeze(mean(mean(Rm < Rth,1),2))';
    oNom(id,:) = ob(1:RB); oAwa(id,:) = ob(RB+1:end);

    %protocol A: best of the first five restarts, scored on Rr -- exp04b's rule
    bn = -inf; ba = -inf; ia = 1; ib = 1;
    for r = 1:5
        en = inf; ea = inf;
        for k = 1:K
            en = min(en, real(THn(:,r)'*Rr(:,:,k)*THn(:,r))/N^2);
            ea = min(ea, real(THa(:,r)'*Rr(:,:,k)*THa(:,r))/N^2);
        end
        if en > bn, bn = en; ia = r; end
        if ea > ba, ba = ea; ib = r; end
    end
    pA(id,:) = [rNom(id,ia) rAwa(id,ib) oNom(id,ia) oAwa(id,ib)];

    fprintf('  D_fb %2d | B: nom %.4f awa %.4f d %+.4f | A: nom %.4f awa %.4f | %.0f s\n', ...
        Dfb, mean(rNom(id,:)), mean(rAwa(id,:)), mean(dRate(id,:)), ...
        pA(id,1), pA(id,2), toc(tA));
    save('exp12b_safe.mat','DfbG','dRate','rNom','rAwa','oNom','oAwa','pA', ...
         'rmsePos','zPred','zSweep','scale','RB','nEval','Rth','rho','wNF','id');
end

%% ---- report -------------------------------------------------------------
se = @(x) 1.96*std(x,0,2)/sqrt(RB);
sd = se(dRate); md = mean(dRate,2); td = md./(sd/1.96);
dO = oAwa - oNom; mo = mean(dO,2); so = se(dO);

fprintf('\n PROTOCOL B (restart-safe): paired min-rate gain, mean +- 95%% CI over %d restarts\n', RB);
fprintf(' D_fb |  nominal   aware   |    gain   +- CI      t   | outage nom  awa\n');
for id = 1:nD
    fprintf(' %4d | %7.4f %7.4f | %+7.4f %7.4f %6.2f | %9.3e %9.3e\n', ...
        DfbG(id), mean(rNom(id,:)), mean(rAwa(id,:)), md(id), sd(id), td(id), ...
        mean(oNom(id,:)), mean(oAwa(id,:)));
end
fprintf('\n resolved above zero at %d of %d operating points\n', sum(md-sd > 0), nD);

fprintf('\n PROTOCOL A (exp04b''s rule, best of 5 scored on Rr) -- reproduction check\n');
for id = 1:nD
    fprintf(' %4d | %7.4f %7.4f | %+7.4f\n', DfbG(id), pA(id,1), pA(id,2), pA(id,2)-pA(id,1));
end

pub = load('exp04b_feedback_sparsity.mat');
fprintf('\n PUBLISHED (exp04b_feedback_sparsity.mat)\n');
for id = 1:nD
    fprintf(' %4d | %7.4f %7.4f | %+7.4f   (dPair %+7.4f)\n', ...
        pub.DfbG(id), pub.minR(id,1), pub.minR(id,2), ...
        pub.minR(id,2)-pub.minR(id,1), pub.dPair(id,1));
end

done = true;
save('exp12b_safe.mat','DfbG','dRate','rNom','rAwa','oNom','oAwa','pA', ...
     'rmsePos','zPred','zSweep','scale','RB','nEval','Rth','rho','wNF','done');
fprintf('\n===== saved exp12b_safe.mat =====\n');

%% ================= local functions =================
function th = localMaxMin(R, th0, pOpt, N, K)
  lo = sqrt(N/K); hi = N; th = functionRobustEqualGainMM(R, th0, lo, pOpt);
  for it = 1:9
    gb = 0.5*(lo+hi);
    t  = functionRobustEqualGainMM(R, th0, gb, pOpt);
    psi = zeros(K,1);
    for k = 1:K, psi(k) = sqrt(max(real(t'*R(:,:,k)*t),0)); end
    if min(psi) >= 0.98*gb, lo = gb; th = t; else, hi = gb; end
  end
end

function Sig = localCov(Ts,d,sa,r)
%Steady-state Kalman prediction covariance, copied verbatim from exp04b so the
%operating points are identical.
F=[eye(3) Ts*eye(3); zeros(3) eye(3)]; H=[eye(3) zeros(3)]; R=r*eye(3);
Q=sa^2*[Ts^4/4*eye(3) Ts^3/2*eye(3); Ts^3/2*eye(3) Ts^2*eye(3)];
P=eye(6);
for i=1:600
  Pm=F*P*F.'+Q; Kg=Pm*H.'/(H*Pm*H.'+R); Pn=(eye(6)-Kg*H)*Pm;
  if norm(Pn-P,'fro')<1e-15, P=Pn; break; end
  P=Pn;
end
for i=1:d, P=F*P*F.'+Q; end
Sig=P(1:3,1:3);
end
