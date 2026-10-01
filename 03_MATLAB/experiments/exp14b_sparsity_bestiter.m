addPaper2Paths();
%exp14b_sparsity_bestiter
%
%SECTION V-F (Fig. 8) ON THE BEST-ITERATE BISECTION.
%
%Supersedes exp12b/exp12c. Everything about the protocol is carried over
%unchanged -- restart bank drawn once, nominal designs solved ONCE and reused at
%every feedback interval, nothing selected, differences formed within a restart
%and within an evaluation draw -- with two changes:
%
%1. THE BISECTION STOPPING RULE. exp12b returns the last iterate ACCEPTED by
%   min_k psi_k >= 0.98*gbar. The equal-gain MM drives psi toward gbar, so
%   lowering gbar lowers the achieved psi with it and good probes are discarded.
%   This script keeps every probe and returns the iterate with the largest
%   min_k psi_k. Strictly at least as good as the published rule, which is
%   computed alongside on the same restarts and the same draws so its cost is
%   measured rather than asserted.
%
%2. RESTART BANK 12 -> 24, matching exp14a so Sections V-E and V-F are quoted
%   at the same precision.
%
%GEOMETRY AND BUDGET. rng(23) with exp04b's draw order, so phat, rho and w_NF
%are identical to exp04b and exp12b. The link-budget scale is LOADED from
%exp12b_safe.mat rather than recomputed, so every absolute rate here is on the
%same budget as the published ones.
%
%RUNTIME. About 84 dense bisections at roughly 5 s each, plus the swept-moment
%assembly, so of order 20 minutes. Do NOT issue MATLAB commands while this runs
%-- they share the command queue and will take the script down. Poll
%exp14b_bestiter.mat from the filesystem instead; it is saved after every row.
%
%Version 1.0 (2026-09-19). License: GPLv2.

clear; close all; rng(23,'twister');
fprintf('\n===== exp14b: V-F on the best-iterate bisection =====\n');

p = functionSimParams(); P_ris = functionRISPositions(p);
N = p.N; K = p.K; lam = p.lambda;
pOpt = p; pOpt.Imax = 150; pOpt.muTol = 1e-4;

phat = zeros(3,K);
for k = 1:K
    phat(:,k) = [p.region(1,1)+diff(p.region(1,:))*rand;
                 p.region(2,1)+diff(p.region(2,:))*rand;
                 p.region(3,1)+diff(p.region(3,:))*rand];
end
rho = mean(vecnorm(phat,2,1));  wNF = lam*rho/p.D;

Ts = 1e-3;  v = 60/3.6;  vhat = [0;1;0];
DfbG = [1 2 5 10 20 50];  nD = numel(DfbG);
nEval = 3000;  RB = 24;  Rth = 2;
TH0 = exp(1j*2*pi*rand(N,RB));
S12 = load('exp12b_safe.mat');  scale = S12.scale;
fprintf('  N=%d K=%d Ts=%.1f ms v=%.1f km/h w_NF=%.5f m | %d restarts, %d draws, scale from exp12b\n', ...
        N, K, Ts*1e3, v*3.6, wNF, RB, nEval);

%% ---- nominal designs: solved ONCE ---------------------------------------
Rn = zeros(N,N,K);
for k = 1:K, Rn(:,:,k) = functionRobustMoment(P_ris, phat(:,k), zeros(3), lam); end
THn = zeros(N,RB);  THnPub = zeros(N,RB);  tA = tic;
for r = 1:RB, [THn(:,r), THnPub(:,r)] = localBisectDense(Rn, TH0(:,r), pOpt, N, K); end
fprintf('  nominal bank solved once (%.0f s)\n', toc(tA));

%% ---- sweep ---------------------------------------------------------------
rNom = zeros(nD,RB); rAwa = zeros(nD,RB); rNomP = zeros(nD,RB); rAwaP = zeros(nD,RB);
oNom = zeros(nD,RB); oAwa = zeros(nD,RB); oNomP = zeros(nD,RB); oAwaP = zeros(nD,RB);
dRate = zeros(nD,RB); dRateP = zeros(nD,RB); dOut = zeros(nD,RB);
rmsePos = zeros(nD,1); zPred = zeros(nD,1); zSweep = zeros(nD,1);

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

    THa = zeros(N,RB);  THaPub = zeros(N,RB);
    for r = 1:RB, [THa(:,r), THaPub(:,r)] = localBisectDense(Rr, TH0(:,r), pOpt, N, K); end

    TH = [THn THa THnPub THaPub];
    es = RandStream('twister','Seed',4242);          % identical draws at every D_fb
    G  = zeros(nEval,K,4*RB);
    for m = 1:nEval
        e  = sigp*randn(es,3,K) + vhat*((rand(es,1,K)-0.5)*delta);
        At = functionArrayResponse(P_ris, phat + e, lam);
        G(m,:,:) = abs(At'*TH).^2/N^2;
    end

    Rm = log2(1 + scale*G);  mn = squeeze(min(Rm,[],2));
    b = @(j) (j-1)*RB + (1:RB);
    rNom(id,:)=mean(mn(:,b(1)),1);  rAwa(id,:)=mean(mn(:,b(2)),1);
    rNomP(id,:)=mean(mn(:,b(3)),1); rAwaP(id,:)=mean(mn(:,b(4)),1);
    dRate(id,:) =mean(mn(:,b(2))-mn(:,b(1)),1);
    dRateP(id,:)=mean(mn(:,b(4))-mn(:,b(3)),1);
    ob = squeeze(mean(mean(Rm < Rth,1),2))';
    oNom(id,:)=ob(b(1)); oAwa(id,:)=ob(b(2)); oNomP(id,:)=ob(b(3)); oAwaP(id,:)=ob(b(4));
    dOut(id,:)=oAwa(id,:)-oNom(id,:);

    fprintf('  D_fb %2d | nom %.4f awa %.4f | gain %+.4f | pub-rule gain %+.4f | %.0f s\n', ...
        Dfb, mean(rNom(id,:)), mean(rAwa(id,:)), mean(dRate(id,:)), mean(dRateP(id,:)), toc(tA));
    save('exp14b_bestiter.mat','DfbG','rNom','rAwa','rNomP','rAwaP','oNom','oAwa', ...
        'oNomP','oAwaP','dRate','dRateP','dOut','rmsePos','zPred','zSweep', ...
        'scale','RB','nEval','Rth','rho','wNF','phat','K','id');
end

%% ---- report --------------------------------------------------------------
sef = @(x) 1.96*std(x,0,2)/sqrt(RB);
tst = @(x) mean(x,2)./(std(x,0,2)/sqrt(RB));

fprintf('\n PAIRED MIN-RATE GAIN, best-iterate bisection, %d restarts\n', RB);
fprintf(' D_fb |  nominal   aware  |    gain    +-CI      t   | %%gain\n');
for i = 1:nD
  fprintf(' %4d | %7.4f %7.4f | %+7.4f %7.4f %6.2f | %+6.1f%%\n', DfbG(i), ...
    mean(rNom(i,:)), mean(rAwa(i,:)), mean(dRate(i,:)), sef(dRate(i,:)), tst(dRate(i,:)), ...
    100*mean(dRate(i,:))/mean(rNom(i,:)));
end
fprintf(' resolved above zero at %d of %d\n', sum(mean(dRate,2)-sef(dRate)>0), nD);

fprintf('\n THE TWO BISECTION RULES, same restarts and draws\n');
fprintf(' D_fb | as-published gain   t   | best-iterate gain   t   | rate lift nom / awa\n');
for i = 1:nD
  fprintf(' %4d | %+9.4f %6.2f | %+9.4f %6.2f | %+6.2f%% %+6.2f%%\n', DfbG(i), ...
    mean(dRateP(i,:)), tst(dRateP(i,:)), mean(dRate(i,:)), tst(dRate(i,:)), ...
    100*(mean(rNom(i,:))-mean(rNomP(i,:)))/mean(rNomP(i,:)), ...
    100*(mean(rAwa(i,:))-mean(rAwaP(i,:)))/mean(rAwaP(i,:)));
end

nS = nEval*K;
fprintf('\n OUTAGE EVENTS per %d user-slots\n', nS);
for i = 1:nD
  fprintf(' %4d | nom %9.1f  awa %9.1f | paired %+9.1f  t %6.2f\n', DfbG(i), ...
    mean(oNom(i,:))*nS, mean(oAwa(i,:))*nS, mean(dOut(i,:))*nS, tst(dOut(i,:)));
end
fprintf('\n restart CV%% nom: %s\n', sprintf('%5.2f ',100*std(rNom,0,2)./mean(rNom,2)));
fprintf(' restart CV%% awa: %s\n', sprintf('%5.2f ',100*std(rAwa,0,2)./mean(rAwa,2)));
fprintf(' RMSE [mm]: %s\n', sprintf('%5.1f ',1e3*rmsePos));
fprintf(' zeta pred: %s\n', sprintf('%6.4f ',zPred));
fprintf(' delta/wNF: %s\n', sprintf('%6.4f ',zSweep));

done = true;
save('exp14b_bestiter.mat','DfbG','rNom','rAwa','rNomP','rAwaP','oNom','oAwa', ...
    'oNomP','oAwaP','dRate','dRateP','dOut','rmsePos','zPred','zSweep', ...
    'scale','RB','nEval','Rth','rho','wNF','phat','K','done');
fprintf('\n===== saved exp14b_bestiter.mat =====\n');

%% ================= local functions =================
function [thBest, thPub] = localBisectDense(R, th0, pOpt, N, K)
  lo = sqrt(N/K); hi = N;
  t0 = functionRobustEqualGainMM(R, th0, lo, pOpt);
  thPub = t0; thBest = t0; best = localMinPsiDense(R, t0, K);
  for it = 1:9
    gb = 0.5*(lo+hi);
    t  = functionRobustEqualGainMM(R, th0, gb, pOpt);
    mp = localMinPsiDense(R, t, K);
    if mp > best, best = mp; thBest = t; end
    if mp >= 0.98*gb, lo = gb; thPub = t; else, hi = gb; end
  end
end
function mp = localMinPsiDense(R, t, K)
  ps = zeros(K,1);
  for k = 1:K, ps(k) = sqrt(max(real(t'*R(:,:,k)*t),0)); end
  mp = min(ps);
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
