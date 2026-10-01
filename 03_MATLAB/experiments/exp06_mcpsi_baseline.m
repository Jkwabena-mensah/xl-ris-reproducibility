%exp06_mcpsi_baseline
%
%PANEL REVIEW ITEM M1.
%
%Proposition 2 rules out sample-averaging the MINIMUM-DISTANCE objective (P0)
%with the auxiliary phase omega frozen across draws. It says nothing about a
%Monte-Carlo estimate of
%
%       psi_k^2(theta) = theta^H Rbar_k theta,
%
%which is invariant to a phase common to all elements BY CONSTRUCTION, and so
%is untouched by the common-mode decoherence that kills (P0). The estimator
%
%       Rhat_k = (1/S) sum_s a(p_s) a(p_s)^H = A_k A_k^H,
%       A_k    = [a(p_1) ... a(p_S)] / sqrt(S)        (N x S)
%
%is itself a FACTOR. It therefore drops straight into functionRobustMaxMinMM's
%existing factored path -- same solver, same warm start, same bisection, same
%softmin annealing -- with NO thin QR and NO eigendecomposition. Only the
%factor differs: rank r = 2 (closed form, compressed) against rank S (MC).
%
%WHAT THIS SETTLES. The paper's negative result is that uncertainty-aware
%design does not pay because the moment construction lands in t_fix, and t_fix
%is what breaks feasibility at short intervals. MC-psi pays no such t_fix. So:
%
%   - if MC-psi is competitive at an affordable S, the negative result is about
%     the FACTORIZATION CHOICE, not about uncertainty-awareness, and the paper
%     must say so;
%   - if MC-psi needs an S so large that O(S K N) per iteration destroys M, or
%     if at affordable S it overfits its own draws, the negative result stands
%     exactly as written and is now defended against the obvious objection.
%
%Either outcome is publishable. Not knowing which is not.
%
%PROTOCOL.
%  (a) Every pipeline is timed in ONE interleaved round -- the exp05b lesson.
%      Wall-clock on a loaded desktop drifts over minutes, so pipelines timed
%      in separate sessions cannot be compared.
%  (b) Every design is scored on INDEPENDENT draws, never on the draws it was
%      designed from (Remark 3). The gap between in-sample and out-of-sample
%      score is the overfitting cost of sampling, and it is the quantity the
%      closed form does not pay.
%  (c) Designs are compared on COMMON position draws and a common warm start,
%      so differences are paired.
%
%This is version 1.0 (Last edited: 2026-09-13)
%License: GPLv2.

clear; close all; rng(606,'twister');
addPaper2Paths();
fprintf('\n================ exp06: Monte-Carlo psi baseline (M1) ================\n');

p     = functionSimParams();
P_ris = functionRISPositions(p);
N = p.N; K = p.K; lam = p.lambda;

% ---- operating point: the one Table IV and Section V-A report
SigR  = (0.01)^2*eye(3);        % 1 cm position uncertainty
delta = 0.05;                   % 5 cm traverse over the hold
vhat  = [0;1;0];                % transverse -- the binding direction
Ts    = 2e-3;                   % reconfiguration interval [s]

pU = zeros(3,K);
for k = 1:K
    pU(:,k) = [p.region(1,1)+diff(p.region(1,:))*rand;
               p.region(2,1)+diff(p.region(2,:))*rand;
               p.region(3,1)+diff(p.region(3,:))*rand];
end
th0 = exp(1j*2*pi*rand(N,1));   % common warm start for every pipeline

Slist = [2 4 8 16 32 64 128 256];
nRepT = 60;                     % interleaved timing rounds
nRepD = 12;                     % independent design replications
nEval = 4000;                   % independent draws for scoring

sq = @(x) x(:);
drawPos = @(k) pU(:,k) + chol(SigR,'lower')*randn(3,1) + (rand-0.5)*delta*vhat;

%% ---------- 1. closed-form reference ----------
fprintf('\n[1] closed-form factor (Thm 2/3, degree 2, J=6, compressed to r=2)\n');
Fcf = cell(1,K); rcf = 0;
for k = 1:K
    [Fcf{k}, inf1] = functionSweptMoment(P_ris,pU(:,k),SigR,delta,vhat,lam,2,2,6);
    rcf = max(rcf, size(Fcf{k},2));
end
fprintf('    working rank r = %d, full Khatri-Rao rank = %d, energy kept %.6f\n', ...
        rcf, inf1.rankFull, inf1.energyKept);

%% ---------- 2. interleaved timing, all pipelines in one round ----------
fprintf('\n[2] interleaved timing, %d rounds (assembly and per-iteration)\n', nRepT);
nS = numel(Slist);
tAsmCF = zeros(nRepT,1); tItCF = zeros(nRepT,1);
tAsmMC = zeros(nRepT,nS); tItMC = zeros(nRepT,nS);

Fmc = cell(1,nS);
for is = 1:nS
    Fmc{is} = cell(1,K);
    for k = 1:K, Fmc{is}{k} = mcFactor(P_ris,pU(:,k),SigR,delta,vhat,lam,Slist(is)); end
end

for r = 1:nRepT
    tic; for k=1:K, Fcf{k} = functionSweptMoment(P_ris,pU(:,k),SigR,delta,vhat,lam,2,2,6); end
    tAsmCF(r) = toc;
    tic; for k=1:K, v = Fcf{k}*(Fcf{k}'*th0); end, tItCF(r) = toc; %#ok<NASGU>
    for is = 1:nS
        S = Slist(is);
        tic; for k=1:K, Fmc{is}{k} = mcFactor(P_ris,pU(:,k),SigR,delta,vhat,lam,S); end
        tAsmMC(r,is) = toc;
        tic; for k=1:K, v = Fmc{is}{k}*(Fmc{is}{k}'*th0); end, tItMC(r,is) = toc; %#ok<NASGU>
    end
end

t_asm_cf = median(tAsmCF);  t_it_cf = median(tItCF);
t_asm_mc = median(tAsmMC,1); t_it_mc = median(tItMC,1);

%% ---------- 3. design and out-of-sample score ----------
fprintf('\n[3] design + independent scoring, %d replications, %d eval draws\n', nRepD, nEval);
pf = p; pf.Imax = 60; pf.muTol = 1e-2;

gCF = zeros(nRepD,1);
gMC = zeros(nRepD,nS); gMCin = zeros(nRepD,nS);
for rep = 1:nRepD
    thCF = functionRobustMaxMinMM(Fcf, th0, pf);
    gCF(rep) = scoreDesign(P_ris,pU,SigR,delta,vhat,lam,thCF,nEval,9000+rep);
    for is = 1:nS
        Fs = cell(1,K);
        for k=1:K, Fs{k} = mcFactor(P_ris,pU(:,k),SigR,delta,vhat,lam,Slist(is)); end
        th = functionRobustMaxMinMM(Fs, th0, pf);
        gMC(rep,is)   = scoreDesign(P_ris,pU,SigR,delta,vhat,lam,th,nEval,9000+rep);
        ik = zeros(K,1);
        for k=1:K, ik(k) = real(th'*(Fs{k}*(Fs{k}'*th)))/N^2; end
        gMCin(rep,is) = min(ik);            % in-sample, on its own S draws
    end
end

%% ---------- 4. report ----------
ci = @(x) 1.96*std(x)/sqrt(numel(x));
fprintf('\n===================== RESULTS =====================\n');
fprintf('operating point: sigma_p = 1 cm, delta = 5 cm transverse, Ts = %.0f ms, K = %d, N = %d\n', ...
        Ts*1e3, K, N);
fprintf('\nCLOSED FORM (r = %d)\n', rcf);
fprintf('  assembly   %8.3f ms    per-iteration %9.4f ms\n', t_asm_cf*1e3, t_it_cf*1e3);
fprintf('  budget M at Ts = %.0f ms : %8.1f\n', Ts*1e3, (Ts-t_asm_cf)/t_it_cf);
fprintf('  out-of-sample min gain   %.6f  +- %.6f\n', mean(gCF), ci(gCF));

fprintf('\nMONTE-CARLO psi\n');
fprintf('    S   asm[ms]  t_it[ms]       M     out-of-sample      in-sample   optimism   vs closed form\n');
for is = 1:nS
    M = (Ts - t_asm_mc(is))/t_it_mc(is);
    fprintf('  %4d  %8.3f %9.4f %8.1f   %.6f+-%.6f  %.6f   %+7.1f%%   %+7.2f%%\n', ...
        Slist(is), t_asm_mc(is)*1e3, t_it_mc(is)*1e3, M, ...
        mean(gMC(:,is)), ci(gMC(:,is)), mean(gMCin(:,is)), ...
        100*(mean(gMCin(:,is))/mean(gMC(:,is))-1), ...
        100*(mean(gMC(:,is))/mean(gCF)-1));
end

fprintf('\nINTERPRETATION GUIDE\n');
fprintf('  "optimism" = in-sample minus out-of-sample, as a %% -- the overfitting the\n');
fprintf('  closed form does not pay. "vs closed form" > 0 means MC-psi wins.\n');
fprintf('  A pipeline is only viable if M > 0 AND it matches the closed form out of sample.\n');

Smatch = NaN;
for is = 1:nS
    d = gMC(:,is) - gCF;                       % paired
    if mean(d) + 1.96*std(d)/sqrt(nRepD) >= 0  % not significantly worse
        Smatch = Slist(is); break
    end
end
if isnan(Smatch)
    fprintf('\n  VERDICT: no S in the sweep matches the closed form out of sample.\n');
else
    isx = find(Slist==Smatch);
    Mm = (Ts - t_asm_mc(isx))/t_it_mc(isx);
    fprintf('\n  VERDICT: smallest matching S = %d, at which M = %.1f (closed form M = %.1f).\n', ...
            Smatch, Mm, (Ts-t_asm_cf)/t_it_cf);
    if Mm <= 0
        fprintf('           M <= 0: MC-psi is infeasible at the S it needs. Negative result STANDS.\n');
    else
        fprintf('           M > 0: MC-psi is feasible at the S it needs. REPORT THIS in the paper.\n');
    end
end

save('exp06_mcpsi.mat','Slist','t_asm_cf','t_it_cf','t_asm_mc','t_it_mc', ...
     'gCF','gMC','gMCin','Ts','nRepT','nRepD','nEval','rcf','p','SigR','delta');
fprintf('\n================ done. saved exp06_mcpsi.mat ================\n');

%% ---------- helpers ----------
function F = mcFactor(P_ris, pk, Sig, delta, vhat, lam, S)
%N x S Monte-Carlo factor. Rhat = F F^H. No QR, no eigendecomposition.
    L = chol(Sig,'lower');
    E = L*randn(3,S) + (rand(1,S)-0.5).*delta.*vhat(:);
    F = functionArrayResponse(P_ris, pk + E, lam) / sqrt(S);
end

function g = scoreDesign(P_ris, pU, Sig, delta, vhat, lam, theta, nEval, seed)
%Out-of-sample min_k E[|a^H theta|^2]/N^2 on INDEPENDENT draws.
    rs = RandStream('twister','Seed',seed);
    K = size(pU,2); N = numel(theta); L = chol(Sig,'lower');
    gk = zeros(K,1);
    for k = 1:K
        E  = L*randn(rs,3,nEval) + (rand(rs,1,nEval)-0.5).*delta.*vhat(:);
        A  = functionArrayResponse(P_ris, pU(:,k) + E, lam);
        gk(k) = mean(abs(A'*theta).^2)/N^2;
    end
    g = min(gk);
end
