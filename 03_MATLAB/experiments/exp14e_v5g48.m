addPaper2Paths();
%exp14e_v5g48
%
%SECTION V-H AT 48 RESTARTS.
%
%WHY. At 24 restarts (exp14c) the line-of-sight CONTROL in Section V-H is
%resolved only at the top of the zeta range, t = 1.94 at zeta = 0.10 and 2.21 at
%0.12, against t = 4.6 to 8.0 under the last-accepted stopping rule on the same
%bank. The control is what licenses the claim that the cascaded-Rician
%non-result is caused by channel randomness rather than by the comparison being
%empty. If the control cannot be resolved, that argument rests on the spread
%ratio alone. Doubling the bank either firms it up or establishes that the
%line-of-sight gain is simply small relative to solver variability at this
%geometry, and either answer is reportable.
%
%THE BANK IS A STRICT SUPERSET OF exp14c's. rng(11) in exp04c's draw order gives
%phat and the first five restarts; rng(777) then fills the rest. MATLAB fills
%rand(N,43) column by column, so its first 19 columns are bit-identical to
%rand(N,19). Restarts 1 to 24 of this run are therefore the same 24 objects
%exp14c used, and the first-24 subset is reported alongside the full bank as a
%reproduction check. If the subset does not reproduce exp14c, something else
%moved and the 48-restart numbers should not be trusted.
%
%Everything else is copied verbatim from exp14c: the evaluation stream
%(Seed 9001), the scatterer construction, the Rician mixture, the per-zeta link
%budget, and both stopping rules carried through the same evaluation.
%
%RUNTIME of order an hour. The file is saved after every zeta, so a partial run
%is usable. Do NOT issue MATLAB commands while this runs.
%
%Version 1.0 (2026-09-19). License: GPLv2.

clear; close all;
fprintf('\n===== exp14e: Section V-H at 48 restarts =====\n');

rng(11,'twister');                      % exp04c's seed, consumed in its order
p = functionSimParams(); P = functionRISPositions(p);
N = p.N; K = p.K; lam = p.lambda;
phat = zeros(3,K);
for k = 1:K
    phat(:,k) = [p.region(1,1)+diff(p.region(1,:))*rand;
                 p.region(2,1)+diff(p.region(2,:))*rand;
                 p.region(3,1)+diff(p.region(3,:))*rand];
end
rho = mean(vecnorm(phat,2,1)); wNF = lam*rho/p.D;
nBank = 5;
TH0 = exp(1j*2*pi*rand(N,nBank));

pOpt = p; pOpt.Imax = 150; pOpt.muTol = 1e-4;
bin = functionArrayResponse(P, p.pBS, lam);
zetaG = [0.02 0.04 0.06 0.08 0.10 0.12]; nZ = numel(zetaG);
nEval = 3000; SNRop = 10; Lfix = p.Lscat; kap = p.kappa;
RB = 48; RB0 = 24;                      % full bank, and exp14c's subset
rng(777,'twister'); THB = [TH0, exp(1j*2*pi*rand(N,RB-nBank))];
fprintf('  N=%d K=%d rho=%.2f w_NF=%.5f | %d restarts (subset %d), %d draws\n', ...
        N, K, rho, wNF, RB, RB0, nEval);

%% ---- nominal designs, both rules -----------------------------------------
Rn = zeros(N,N,K);
for k = 1:K, Rn(:,:,k) = functionRobustMoment(P, phat(:,k), zeros(3), lam); end
TN = zeros(N,RB); TNp = zeros(N,RB); tA = tic;
for r = 1:RB, [TN(:,r), TNp(:,r)] = localBisectDense(Rn, THB(:,r), pOpt, N, K); end
fprintf('  nominal bank solved (%.0f s)\n', toc(tA));

VL = zeros(nZ,RB); VR = zeros(nZ,RB); VLp = zeros(nZ,RB); VRp = zeros(nZ,RB);
mL = zeros(nZ,1); mR = zeros(nZ,1); tL = zeros(nZ,1); tR = zeros(nZ,1);
mLp = zeros(nZ,1); mRp = zeros(nZ,1); tLp = zeros(nZ,1); tRp = zeros(nZ,1);
tL0 = zeros(nZ,1); tR0 = zeros(nZ,1);   % first-24 subset, reproduction check

for iz = 1:nZ
  sig = zetaG(iz)*wNF; Sig = sig^2*eye(3);
  Rr = zeros(N,N,K);
  for k = 1:K, Rr(:,:,k) = functionRobustMoment(P, phat(:,k), Sig, lam); end
  TU = zeros(N,RB); TUp = zeros(N,RB);
  for r = 1:RB, [TU(:,r), TUp(:,r)] = localBisectDense(Rr, THB(:,r), pOpt, N, K); end

  TT = [TN TU TNp TUp];                        % N x 4RB
  es = RandStream('twister','Seed',9001);
  gL = zeros(nEval,K,4*RB); gR = zeros(nEval,K,4*RB);
  for m = 1:nEval
    Pt = phat + sig*randn(es,3,K);
    for k = 1:K
      Q  = [Pt(:,k), Pt(:,k)+p.scatSpread*[randn(es,2,Lfix); 0.3*randn(es,1,Lfix)]];
      Aa = functionArrayResponse(P,Q,lam); aL = Aa(:,1);
      be = (randn(es,Lfix,1)+1i*randn(es,Lfix,1))/sqrt(2);
      bN = Aa(:,2:end)*be/sqrt(Lfix);
      gL(m,k,:) = abs(aL'*TT).^2/N^2;
      h = conj(sqrt(kap/(kap+1))*aL + sqrt(1/(kap+1))*bN).*bin;
      gR(m,k,:) = abs(h'*TT).^2/N^2;
    end
  end

  VL(iz,:)  = localPairs(gL, RB, 0, SNRop);
  VR(iz,:)  = localPairs(gR, RB, 0, SNRop);
  VLp(iz,:) = localPairs(gL, RB, 2*RB, SNRop);
  VRp(iz,:) = localPairs(gR, RB, 2*RB, SNRop);
  mL(iz)=mean(VL(iz,:));   tL(iz)=mL(iz)/(std(VL(iz,:))/sqrt(RB));
  mR(iz)=mean(VR(iz,:));   tR(iz)=mR(iz)/(std(VR(iz,:))/sqrt(RB));
  mLp(iz)=mean(VLp(iz,:)); tLp(iz)=mLp(iz)/(std(VLp(iz,:))/sqrt(RB));
  mRp(iz)=mean(VRp(iz,:)); tRp(iz)=mRp(iz)/(std(VRp(iz,:))/sqrt(RB));
  v0=VL(iz,1:RB0); tL0(iz)=mean(v0)/(std(v0)/sqrt(RB0));
  w0=VR(iz,1:RB0); tR0(iz)=mean(w0)/(std(w0)/sqrt(RB0));

  fprintf('  zeta %.2f | LoS %+0.4f (t %5.2f) Ric %+0.4f (t %5.2f) | subset-24 t: LoS %5.2f Ric %5.2f | %.0f s\n', ...
     zetaG(iz), mL(iz),tL(iz), mR(iz),tR(iz), tL0(iz),tR0(iz), toc(tA));
  save('exp14e_v5g48.mat','zetaG','VL','VR','VLp','VRp','mL','mR','tL','tR', ...
       'mLp','mRp','tLp','tRp','tL0','tR0','RB','RB0','nEval','rho','wNF','phat','kap','Lfix','iz');
end

%% ---- report --------------------------------------------------------------
sef = @(v) 1.96*std(v)/sqrt(RB);
fprintf('\n PAIRED MIN-RATE GAIN, best-iterate bisection, %d restarts\n', RB);
fprintf(' zeta |  line-of-sight   +-CI      t   |   cascaded Rician  +-CI      t\n');
for i = 1:nZ
  fprintf(' %.2f | %+14.4f %7.4f %6.2f | %+16.4f %7.4f %6.2f\n', zetaG(i), ...
    mL(i), sef(VL(i,:)), tL(i), mR(i), sef(VR(i,:)), tR(i));
end
fprintf(' LoS resolved at %d of %d | Rician resolved at %d of %d | max |t| Rician = %.2f\n', ...
  sum(abs(tL)>1.96), nZ, sum(abs(tR)>1.96), nZ, max(abs(tR)));

fprintf('\n REPRODUCTION CHECK: first %d restarts of this bank vs the full %d\n', RB0, RB);
fprintf(' zeta | LoS t (24) | LoS t (48) | Ric t (24) | Ric t (48)\n');
for i = 1:nZ
  fprintf(' %.2f | %10.2f | %10.2f | %10.2f | %10.2f\n', zetaG(i), tL0(i), tL(i), tR0(i), tR(i));
end

fprintf('\n THE TWO BISECTION RULES at %d restarts\n', RB);
fprintf(' zeta | pub LoS (t)      | best LoS (t)     | pub Ric (t)      | best Ric (t)\n');
for i = 1:nZ
  fprintf(' %.2f | %+8.4f (%5.2f) | %+8.4f (%5.2f) | %+8.4f (%5.2f) | %+8.4f (%5.2f)\n', ...
    zetaG(i), mLp(i),tLp(i), mL(i),tL(i), mRp(i),tRp(i), mR(i),tR(i));
end
fprintf('\n spread ratio (Rician sd / LoS sd), best-iterate: %s\n', sprintf('%5.2f ', std(VR,0,2)./std(VL,0,2)));
fprintf(' spread ratio (Rician sd / LoS sd), as published: %s\n', sprintf('%5.2f ', std(VRp,0,2)./std(VLp,0,2)));
fprintf(' sd inflation from the stopping rule (best / pub), LoS: %s\n', sprintf('%5.2f ', std(VL,0,2)./std(VLp,0,2)));

done = true;
save('exp14e_v5g48.mat','zetaG','VL','VR','VLp','VRp','mL','mR','tL','tR', ...
     'mLp','mRp','tLp','tRp','tL0','tR0','RB','RB0','nEval','rho','wNF','phat','kap','Lfix','done');
fprintf('\n===== saved exp14e_v5g48.mat =====\n');

%% ================= local functions =================
function v = localPairs(g, R, off, SNRop)
  gRef = median(reshape(g(:,:,off+(1:R)),[],1)); c = 10^(SNRop/10)/gRef;
  v = zeros(1,R);
  for r = 1:R
    v(r) = mean( min(log2(1+c*g(:,:,off+R+r)),[],2) - min(log2(1+c*g(:,:,off+r)),[],2) );
  end
end
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
