addPaper2Paths();
%exp12c_finish
%
%Completes exp12b at D_fb = 20 and 50. The original run was INTERRUPTED at
%id = 4, not slowed: a timing probe measured one equal-gain call at D_fb = 20 at
%0.41 s, so twelve restarts there cost under a minute. What stopped it was
%issuing MATLAB commands from outside while the script was running -- those
%execute in the same command queue and take the script down with them. Poll the
%saved .mat from the filesystem instead.
%
%The geometry, the restart bank and the link-budget scale are reproduced exactly:
%rng(23) with exp04b's draw order regenerates phat and TH0, and `scale` is taken
%from the partial exp12b_safe.mat rather than recomputed, so the completed rows
%and the new ones are on one budget.
%
%Version 1.0 (2026-09-19). License: GPLv2.

clear; close all; rng(23,'twister');
fprintf('\n===== exp12c: completing exp12b at D_fb = 20 and 50 =====\n');

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

Ts = 1e-3; v = 60/3.6; vhat = [0;1;0];
DfbG = [1 2 5 10 20 50]; nD = numel(DfbG);
nEval = 3000; RB = 12; Rth = 2;
TH0 = exp(1j*2*pi*rand(N,RB));

S = load('exp12b_safe.mat');
assert(S.id == 4, 'expected the partial file to stop at id = 4');
scale = S.scale;
dRate = S.dRate; rNom = S.rNom; rAwa = S.rAwa; oNom = S.oNom; oAwa = S.oAwa;
rmsePos = S.rmsePos; zPred = S.zPred; zSweep = S.zSweep; pA = S.pA;
fprintf('  resumed: rows 1-4 kept, scale = %.6e\n', scale);

Rn = zeros(N,N,K);
for k = 1:K, Rn(:,:,k) = functionRobustMoment(P_ris, phat(:,k), zeros(3), lam); end
THn = zeros(N,RB);
tA = tic;
for r = 1:RB, THn(:,r) = localMaxMin(Rn, TH0(:,r), pOpt, N, K); end
fprintf('  nominal bank re-solved (%.0f s)\n', toc(tA));

for id = 5:6
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
    es = RandStream('twister','Seed',4242);
    G  = zeros(nEval,K,2*RB);
    for m = 1:nEval
        e  = sigp*randn(es,3,K) + vhat*((rand(es,1,K)-0.5)*delta);
        At = functionArrayResponse(P_ris, phat + e, lam);
        G(m,:,:) = abs(At'*TH).^2/N^2;
    end

    Rm = log2(1 + scale*G);
    mn = squeeze(min(Rm,[],2));
    rNom(id,:)  = mean(mn(:,1:RB),1);
    rAwa(id,:)  = mean(mn(:,RB+1:end),1);
    dRate(id,:) = mean(mn(:,RB+1:end) - mn(:,1:RB), 1);
    ob = squeeze(mean(mean(Rm < Rth,1),2))';
    oNom(id,:) = ob(1:RB); oAwa(id,:) = ob(RB+1:end);

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

    fprintf('  D_fb %2d | nom %.4f awa %.4f | gain %+.4f | A-gain %+.4f | %.0f s\n', ...
        Dfb, mean(rNom(id,:)), mean(rAwa(id,:)), mean(dRate(id,:)), ...
        pA(id,2)-pA(id,1), toc(tA));
    save('exp12b_safe.mat','DfbG','dRate','rNom','rAwa','oNom','oAwa','pA', ...
         'rmsePos','zPred','zSweep','scale','RB','nEval','Rth','rho','wNF','id');
end

%% ---- full report --------------------------------------------------------
sefun = @(x) 1.96*std(x,0,2)/sqrt(RB);
md = mean(dRate,2); sd = sefun(dRate); td = md./(sd/1.96);
pub = load('exp04b_feedback_sparsity_PUBLISHED_backup.mat');

fprintf('\n PROTOCOL B (restart-safe, %d restarts)\n', RB);
fprintf(' D_fb |  nominal   aware  |    gain    +- CI      t   | %%gain\n');
for i = 1:nD
    fprintf(' %4d | %7.4f %7.4f | %+7.4f %7.4f %6.2f | %+6.1f%%\n', ...
        DfbG(i), mean(rNom(i,:)), mean(rAwa(i,:)), md(i), sd(i), td(i), ...
        100*md(i)/mean(rNom(i,:)));
end
fprintf('\n resolved above zero at %d of %d\n', sum(md-sd>0), nD);

fprintf('\n D_fb | A:best-of-5 | published dPair | B:restart-safe | pub %%   | safe %%\n');
for i = 1:nD
    fprintf(' %4d | %+7.4f     | %+7.4f         | %+7.4f        | %+6.1f%% | %+6.1f%%\n', ...
        DfbG(i), pA(i,2)-pA(i,1), pub.dPair(i,1), md(i), ...
        100*(pub.minR(i,2)-pub.minR(i,1))/pub.minR(i,1), 100*md(i)/mean(rNom(i,:)));
end
fprintf('\n restart CV%% nominal: %s\n', sprintf('%5.2f ',100*std(rNom,0,2)./mean(rNom,2)));
fprintf(' restart CV%% aware  : %s\n', sprintf('%5.2f ',100*std(rAwa,0,2)./mean(rAwa,2)));

done = true;
save('exp12b_safe.mat','DfbG','dRate','rNom','rAwa','oNom','oAwa','pA', ...
     'rmsePos','zPred','zSweep','scale','RB','nEval','Rth','rho','wNF','done');
fprintf('\n===== exp12b_safe.mat complete =====\n');

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
