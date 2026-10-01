addPaper2Paths();
%exp14a_rate_outage_bestiter
%
%SECTION V-E (Figs 6 and 7) AND REVIEWER ITEM C-2, ON ONE PROTOCOL.
%
%Supersedes exp12a (Section V-E) and exp13b (the worst-case baseline). Three
%designs are solved on a common 24-restart bank and evaluated on common draws:
%
%  nominal      -- second moment at Sigma = 0
%  aware        -- closed-form second moment at Sigma = (zeta w_NF)^2 I
%  worst-case   -- M = 128 boundary points of the chi^2_{3,0.95} confidence
%                  ellipsoid U_k, giving KM = 1024 rank-one amplitude
%                  constraints in place of K = 8 second moments
%
%TWO CHANGES FROM exp12a, BOTH DELIBERATE.
%
%1. THE BISECTION STOPPING RULE. exp12a's bisection accepts a probe only if
%   min_k psi_k >= 0.98*gbar and returns the last ACCEPTED iterate. The
%   equal-gain MM drives every psi toward gbar, so lowering gbar lowers the
%   achieved psi with it; with clustered constraints the achieved minimum sits a
%   few per cent below gbar and good probes are discarded. This script keeps
%   every probe and returns the iterate with the largest min_k psi_k. It is a
%   STRICT improvement: the returned design is by construction at least as good
%   as the one the published rule returns.
%
%   The published rule is computed ALONGSIDE, on the same restarts and the same
%   draws, so its cost is measured rather than asserted. Nothing else differs
%   between the two columns.
%
%2. RESTART BANK 12 -> 24, and M 64 -> 128. exp13b showed the worst-case design
%   has a wider restart spread (CV up to 7.4%) than the other two and a
%   non-monotone mean-rate difference, which is the signature of the rank-one
%   bisection occasionally settling on a poorer local optimum. Both are raised
%   here before anything is drawn.
%
%WHAT IS UNCHANGED. rng(11) and exp04a's draw order, so phat, rho and w_NF are
%identical to exp12a and exp04a. The link-budget scale is LOADED from
%exp12a_safe.mat rather than recomputed, so every absolute rate in this file is
%on the same budget as the published ones and the comparison is like for like.
%Nominal designs solved ONCE and reused at every zeta; nothing selected
%anywhere; differences formed within a restart and within an evaluation draw.
%
%Version 1.0 (2026-09-19). License: GPLv2.

clear; close all; rng(11,'twister');
fprintf('\n===== exp14a: V-E and C-2 on the best-iterate bisection =====\n');

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

zetaG = [0.02 0.04 0.06 0.08 0.10 0.12];  nZ = numel(zetaG);
nEval = 3000;  RB = 24;  Rth = 2;  chi2 = 7.814728;  M = 128;
TH0 = exp(1j*2*pi*rand(N,RB));
S12 = load('exp12a_safe.mat');  scale = S12.scale;
fprintf('  N=%d K=%d rho=%.2f w_NF=%.5f | %d restarts, %d draws, M=%d, scale from exp12a\n', ...
        N, K, rho, wNF, RB, nEval, M);

%% ---- nominal designs: solved ONCE ---------------------------------------
Rn = zeros(N,N,K);
for k = 1:K, Rn(:,:,k) = functionRobustMoment(P_ris, phat(:,k), zeros(3), lam); end
THn = zeros(N,RB);  THnPub = zeros(N,RB);  tA = tic;
for r = 1:RB, [THn(:,r), THnPub(:,r)] = localBisectDense(Rn, TH0(:,r), pOpt, N, K); end
fprintf('  nominal bank solved once (%.0f s)\n', toc(tA));

%% ---- sweep ---------------------------------------------------------------
U = localFib(M);
rNom  = zeros(nZ,RB); rAwa = zeros(nZ,RB); rWc = zeros(nZ,RB);
rNomP = zeros(nZ,RB); rAwaP = zeros(nZ,RB);
oNom  = zeros(nZ,RB); oAwa = zeros(nZ,RB); oWc = zeros(nZ,RB);
oNomP = zeros(nZ,RB); oAwaP = zeros(nZ,RB);
dAN = zeros(nZ,RB); dWN = zeros(nZ,RB); dAW = zeros(nZ,RB);
dON = zeros(nZ,RB); dOW = zeros(nZ,RB); dANP = zeros(nZ,RB); dONP = zeros(nZ,RB);
tAsm = zeros(nZ,2);

for iz = 1:nZ
    sig = zetaG(iz)*wNF;  Sig = sig^2*eye(3);

    tic; Rr = zeros(N,N,K);
    for k = 1:K, Rr(:,:,k) = functionRobustMoment(P_ris, phat(:,k), Sig, lam); end
    tAsm(iz,1) = toc;

    tic; Lz = chol(Sig,'lower');  Dz = sqrt(chi2)*(Lz*U.');
    Aw = zeros(N,K*M); c = 0;
    for k = 1:K, Aw(:,c+(1:M)) = functionArrayResponse(P_ris, phat(:,k)+Dz, lam); c = c+M; end
    tAsm(iz,2) = toc;

    THa = zeros(N,RB); THaPub = zeros(N,RB); THw = zeros(N,RB);
    for r = 1:RB
        [THa(:,r), THaPub(:,r)] = localBisectDense(Rr, TH0(:,r), pOpt, N, K);
        THw(:,r)                = localBisectRank1(Aw, TH0(:,r), pOpt, N, K);
    end

    TH = [THn THa THw THnPub THaPub];
    es = RandStream('twister','Seed',7000+iz);
    G  = zeros(nEval,K,5*RB);
    for m = 1:nEval
        At = functionArrayResponse(P_ris, phat + sig*randn(es,3,K), lam);
        G(m,:,:) = abs(At'*TH).^2/N^2;
    end
    Rm = log2(1 + scale*G);  mn = squeeze(min(Rm,[],2));
    b = @(j) (j-1)*RB + (1:RB);
    rNom(iz,:)=mean(mn(:,b(1)),1);  rAwa(iz,:)=mean(mn(:,b(2)),1);  rWc(iz,:)=mean(mn(:,b(3)),1);
    rNomP(iz,:)=mean(mn(:,b(4)),1); rAwaP(iz,:)=mean(mn(:,b(5)),1);
    dAN(iz,:)=mean(mn(:,b(2))-mn(:,b(1)),1);
    dWN(iz,:)=mean(mn(:,b(3))-mn(:,b(1)),1);
    dAW(iz,:)=mean(mn(:,b(2))-mn(:,b(3)),1);
    dANP(iz,:)=mean(mn(:,b(5))-mn(:,b(4)),1);

    ob = squeeze(mean(mean(Rm < Rth,1),2))';
    oNom(iz,:)=ob(b(1)); oAwa(iz,:)=ob(b(2)); oWc(iz,:)=ob(b(3));
    oNomP(iz,:)=ob(b(4)); oAwaP(iz,:)=ob(b(5));
    dON(iz,:)=oAwa(iz,:)-oNom(iz,:);  dOW(iz,:)=oWc(iz,:)-oNom(iz,:);
    dONP(iz,:)=oAwaP(iz,:)-oNomP(iz,:);

    fprintf('  zeta %.2f | nom %.4f awa %.4f wc %.4f | a-n %+.4f w-n %+.4f | pub nom %.4f awa %.4f | %.0f s\n', ...
        zetaG(iz), mean(rNom(iz,:)), mean(rAwa(iz,:)), mean(rWc(iz,:)), ...
        mean(dAN(iz,:)), mean(dWN(iz,:)), mean(rNomP(iz,:)), mean(rAwaP(iz,:)), toc(tA));
    save('exp14a_bestiter.mat','zetaG','rNom','rAwa','rWc','rNomP','rAwaP', ...
        'oNom','oAwa','oWc','oNomP','oAwaP','dAN','dWN','dAW','dON','dOW','dANP','dONP', ...
        'tAsm','chi2','M','RB','nEval','Rth','scale','rho','wNF','phat','K','iz');
end

%% ---- report --------------------------------------------------------------
sef = @(x) 1.96*std(x,0,2)/sqrt(RB);
tst = @(x) mean(x,2)./(std(x,0,2)/sqrt(RB));
nS  = nEval*K;

fprintf('\n DELIVERED MIN RATE [bit/s/Hz], best-iterate bisection, %d restarts\n', RB);
fprintf(' zeta | nominal   aware   worst-case |  aware-nom   +-CI      t   |  worst-nom   +-CI      t\n');
for i = 1:nZ
  fprintf(' %.2f | %7.4f %7.4f %9.4f | %+9.4f %7.4f %6.2f | %+9.4f %7.4f %6.2f\n', zetaG(i), ...
    mean(rNom(i,:)), mean(rAwa(i,:)), mean(rWc(i,:)), ...
    mean(dAN(i,:)), sef(dAN(i,:)), tst(dAN(i,:)), ...
    mean(dWN(i,:)), sef(dWN(i,:)), tst(dWN(i,:)));
end
fprintf(' aware-nom resolved above zero at %d of %d\n', sum(mean(dAN,2)-sef(dAN)>0), nZ);

fprintf('\n THE TWO BISECTION RULES, same restarts and draws\n');
fprintf(' zeta | as-published gain  t   | best-iterate gain  t   | rate lift nom / awa\n');
for i = 1:nZ
  fprintf(' %.2f | %+9.4f %6.2f | %+9.4f %6.2f | %+6.2f%% %+6.2f%%\n', zetaG(i), ...
    mean(dANP(i,:)), tst(dANP(i,:)), mean(dAN(i,:)), tst(dAN(i,:)), ...
    100*(mean(rNom(i,:))-mean(rNomP(i,:)))/mean(rNomP(i,:)), ...
    100*(mean(rAwa(i,:))-mean(rAwaP(i,:)))/mean(rAwaP(i,:)));
end

fprintf('\n OUTAGE EVENTS per %d user-slots (threshold %.1f bit/s/Hz)\n', nS, Rth);
fprintf(' zeta |  nominal     aware  worst-case |  awa-nom     t   |  wc-nom      t   | wc vs nom\n');
for i = 1:nZ
  red = NaN; if mean(oNom(i,:))>0, red = 100*(1-mean(oWc(i,:))/mean(oNom(i,:))); end
  fprintf(' %.2f | %9.1f %9.1f %9.1f | %+9.1f %6.2f | %+9.1f %6.2f | %6.1f%%\n', zetaG(i), ...
    mean(oNom(i,:))*nS, mean(oAwa(i,:))*nS, mean(oWc(i,:))*nS, ...
    mean(dON(i,:))*nS, tst(dON(i,:)), mean(dOW(i,:))*nS, tst(dOW(i,:)), red);
end

fprintf('\n ASSEMBLY [ms] closed-form moment vs %d-point ellipsoid\n', M);
for i = 1:nZ
  fprintf(' %.2f | %7.3f %7.3f  ratio %5.2f\n', zetaG(i), 1e3*tAsm(i,1), 1e3*tAsm(i,2), tAsm(i,2)/tAsm(i,1));
end
fprintf('\n restart CV%%  nom: %s\n', sprintf('%5.2f ',100*std(rNom,0,2)./mean(rNom,2)));
fprintf(' restart CV%%  awa: %s\n', sprintf('%5.2f ',100*std(rAwa,0,2)./mean(rAwa,2)));
fprintf(' restart CV%%  wc : %s\n', sprintf('%5.2f ',100*std(rWc ,0,2)./mean(rWc ,2)));

done = true;
save('exp14a_bestiter.mat','zetaG','rNom','rAwa','rWc','rNomP','rAwaP', ...
    'oNom','oAwa','oWc','oNomP','oAwaP','dAN','dWN','dAW','dON','dOW','dANP','dONP', ...
    'tAsm','chi2','M','RB','nEval','Rth','scale','rho','wNF','phat','K','done');
fprintf('\n===== saved exp14a_bestiter.mat =====\n');

%% ================= local functions =================
function [thBest, thPub] = localBisectDense(R, th0, pOpt, N, K)
%Max-min by bisection on the common target (F-35). Returns BOTH the best iterate
%probed and the one the published acceptance rule would have returned.
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
function thBest = localBisectRank1(A, th0, pOpt, N, K)
%Same bisection over the KM rank-one amplitude constraints |a_j' theta|.
  lo = sqrt(N/K); hi = N; Rsum = A*A'; lamBar = max(sum(abs(Rsum),2));
  t0 = localRank1MM(A, th0, lo, pOpt.Imax, Rsum, lamBar);
  thBest = t0; best = min(abs(A'*t0));
  for it = 1:9
    gb = 0.5*(lo+hi);
    t  = localRank1MM(A, th0, gb, pOpt.Imax, Rsum, lamBar);
    mp = min(abs(A'*t));
    if mp > best, best = mp; thBest = t; end
    if mp >= 0.98*gb, lo = gb; else, hi = gb; end
  end
end
function theta = localRank1MM(A, thetaInit, gbar, Imax, Rsum, lamBar)
%Equal-gain MM with rank-one constraints: the Sigma = 0 reduction of
%functionRobustEqualGainMM applied to each ellipsoid boundary point.
  theta = thetaInit(:); theta = theta./abs(theta);
  for t = 1:Imax
    c = A'*theta; psi = abs(c);
    v = lamBar*theta - Rsum*theta + gbar*(A*(c./max(psi,eps)));
    tn = exp(1j*angle(v));
    if norm(tn-theta) < 1e-4 && t > 5, theta = tn; break; end
    theta = tn;
  end
end
function U = localFib(M)
%Fibonacci sphere: M near-uniform directions on the unit sphere.
  i = (0:M-1)'+0.5; phi = acos(1-2*i/M); g = pi*(1+sqrt(5)); t = g*i;
  U = [cos(t).*sin(phi), sin(t).*sin(phi), cos(phi)];
end
