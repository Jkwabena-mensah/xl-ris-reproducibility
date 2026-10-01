addPaper2Paths();
%exp04b_feedback_sparsity
%
%Closes three gaps identified by auditing the manuscript against the PhD-by-
%publication plan for Paper 2, whose central research question asks how the design
%behaves when position estimates are "inaccurate, delayed, correlated, or
%UNAVAILABLE for some time slots".
%
%  (i)   POSITION-PREDICTION RMSE. Declared as a performance measure in
%        Section II-H and never reported. It is reported here as a function of how
%        often position reports actually arrive.
%  (ii)  FEEDBACK SPARSITY / DELAY. Section V asserts that the design begins to pay
%        once feedback becomes sparser than about one report every ten intervals,
%        but never sweeps it. Sparsity is the operationally meaningful axis: it is
%        what a deployment controls, whereas zeta is derived.
%  (iii) COLD START vs WARM START. The plan asks for a comparison against cold-start
%        BCD-MM; Section V reports warm-start iteration counts with nothing to
%        compare them to.
%
%WHY SPARSITY DRIVES BOTH EFFECTS. Dropping reports inflates the prediction
%covariance through the propagation of Appendix A, AND lengthens the displacement
%swept while one configuration is held, delta = ||v|| Ts D_fb. Theorem 2 consumes
%both through the product kernel, so the uncertainty-aware design here is built
%with functionSweptMoment, not with the prediction covariance alone. Travel is
%transverse, which Remark 3 identifies as the binding direction; radial travel is
%nearly free and would understate the effect (this is the F-26 trap).
%
%Both designs are solved by the max-min bisection of Algorithm 1 (see F-35: calling
%the equal-gain iteration at a fixed target is a feasibility test, not an optimiser,
%and reverses the sign of the comparison).
%
%This is version 1.0 (Last edited: 2026-09-12)
%License: GPLv2.

clear; close all; rng(23,'twister');
fprintf('exp04b: feedback sparsity, prediction RMSE, cold vs warm start\n');

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

Ts   = 1e-3;                 % reconfiguration interval [s]
v    = 60/3.6;               % speed [m/s]
vhat = [0;1;0];              % transverse travel, the binding direction
DfbG = [1 2 5 10 20 50];     % position reports every D_fb intervals
nD   = numel(DfbG);
nEval = 3000; nRest = 5; Rth = 2; SNRop = 10;

%COMMON RANDOM NUMBERS (added after the exp04c control exposed the confound).
%The nominal design never sees D_fb, so restart-to-restart variability of the
%nonconvex solver must not be allowed to vary along the sweep. The restart bank
%and the evaluation draws are therefore held identical at every D_fb.
TH0 = exp(1j*2*pi*rand(N,nRest));

fprintf('  N = %d, K = %d, Ts = %.1f ms, v = %.1f km/h, w_NF = %.4f m\n', ...
        N, K, Ts*1e3, v*3.6, wNF);

gammaAll = zeros(nEval,K,nD,2);
rmsePos  = zeros(nD,1); zPred = zeros(nD,1); zSweep = zeros(nD,1);
itCold   = zeros(nD,1); itWarm = zeros(nD,1);

thPrev = exp(1j*2*pi*rand(N,1));
tAll = tic;
for id = 1:nD
    Dfb   = DfbG(id);
    Sig   = localCov(Ts, Dfb, 1.0, p.r);
    sigp  = sqrt(trace(Sig)/3);
    delta = v*Ts*Dfb;
    rmsePos(id) = sqrt(trace(Sig));          % 3-D position-prediction RMSE [m]
    zPred(id)   = sigp/wNF;
    zSweep(id)  = delta/wNF;

    Rn = zeros(N,N,K); Rr = zeros(N,N,K);
    for k = 1:K
        Rn(:,:,k) = functionRobustMoment(P_ris, phat(:,k), zeros(3), lam);
        Fk = functionSweptMoment(P_ris, phat(:,k), Sig, delta, vhat, lam, 4, 2, 6);
        Rr(:,:,k) = Fk*Fk';
    end

    bn = -inf; br = -inf; an = []; ar = []; gbBest = sqrt(N/K);
    for r = 1:nRest
        th0 = TH0(:,r);                     % common restarts across D_fb
        [tn,~]      = localMaxMin(Rn, th0, pOpt, N, K);
        [tr, gbr]   = localMaxMin(Rr, th0, pOpt, N, K);
        en = inf; er = inf;
        for k = 1:K
            en = min(en, real(tn'*Rr(:,:,k)*tn)/N^2);
            er = min(er, real(tr'*Rr(:,:,k)*tr)/N^2);
        end
        if en > bn, bn = en; an = tn; end
        if er > br, br = er; ar = tr; gbBest = gbr; end
    end

    % ---- cold start vs warm start, at the feasible target the bisection found
    [~, itC] = functionRobustEqualGainMM(Rr, exp(1j*2*pi*rand(N,1)), gbBest, pOpt);
    [~, itW] = functionRobustEqualGainMM(Rr, thPrev,                 gbBest, pOpt);
    itCold(id) = itC; itWarm(id) = itW; thPrev = ar;

    evalStream = RandStream('twister','Seed',4242);   % identical draws at every D_fb
    for m = 1:nEval
        e  = sigp*randn(evalStream,3,K) + vhat*((rand(evalStream,1,K)-0.5)*delta);
        At = functionArrayResponse(P_ris, phat + e, lam);
        gammaAll(m,:,id,1) = abs(At'*an).^2.'/N^2;
        gammaAll(m,:,id,2) = abs(At'*ar).^2.'/N^2;
    end
    fprintf('  D_fb %2d | RMSE %6.1f mm | z_pred %.3f | z_sweep %.3f | cold %3d warm %3d | %.0f s\n', ...
        Dfb, 1e3*rmsePos(id), zPred(id), zSweep(id), itC, itW, toc(tAll));
end

%% ---- link budget referenced to the nominal design at D_fb = 1
gRef  = median(reshape(gammaAll(:,:,1,1), [], 1));
scale = 10^(SNRop/10)/gRef;
R     = log2(1 + scale*gammaAll);
minPD = min(R, [], 2);
minR  = squeeze(mean(minPD,1));
minCI = squeeze(1.96*std(minPD,0,1)/sqrt(nEval));
pOut  = squeeze(mean(mean(R < Rth, 2), 1));
outCI = 1.96*sqrt(max(pOut.*(1-pOut),0)/(nEval*K));

fprintf('\n D_fb | RMSE [mm] | z_pred | z_sweep | min rate nom/aware (+-CI) | P_out nom/aware\n');
for id = 1:nD
    fprintf(' %4d | %8.1f  | %.3f  | %.3f   | %5.2f / %5.2f (+-%.3f)      | %.4f / %.4f\n', ...
        DfbG(id), 1e3*rmsePos(id), zPred(id), zSweep(id), ...
        minR(id,1), minR(id,2), mean(minCI(id,:)), pOut(id,1), pOut(id,2));
end
dPair = zeros(nD,3); nEv = zeros(nD,2);
for id = 1:nD
    Rz = log2(1 + scale*squeeze(gammaAll(:,:,id,:)));
    dd = min(Rz(:,:,2),[],2) - min(Rz(:,:,1),[],2);        % paired per draw
    dPair(id,:) = [mean(dd), mean(dd)-1.96*std(dd)/sqrt(nEval), mean(dd)+1.96*std(dd)/sqrt(nEval)];
    nEv(id,:)   = [sum(sum(Rz(:,:,1)<Rth)), sum(sum(Rz(:,:,2)<Rth))];
end
fprintf('\nPAIRED difference (aware - nominal) and outage EVENT COUNTS out of %d user-slots:\n', nEval*K);
for id = 1:nD
    fprintf(' D_fb %2d | %+7.4f [%+.4f,%+.4f] %-5s | outage %5d -> %5d\n', ...
        DfbG(id), dPair(id,1), dPair(id,2), dPair(id,3), string(dPair(id,2)>0), nEv(id,1), nEv(id,2));
end
sig95 = dPair(:,2) > 0;
fprintf('\n significant at 95%%: %d of %d | cold/warm iteration ratio %.2f (mean)\n', ...
    sum(sig95), nD, mean(itCold./max(itWarm,1)));
fprintf(' prediction RMSE grows %.1f mm -> %.1f mm over the sweep; sweep ratio dominates from D_fb = %d\n', ...
    1e3*rmsePos(1), 1e3*rmsePos(end), DfbG(find(zSweep > zPred, 1)));

%% ---- figure
figdir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'06_Figures');
if ~isfolder(figdir), mkdir(figdir); end
fh = figure('Units','centimeters','Position',[2 2 8.8 6.8],'Color','w');
hold on; box on; grid on;
errorbar(DfbG, minR(:,1), minCI(:,1), '-o','Color',[0 0 0], ...
    'LineWidth',1.4,'MarkerSize',5,'CapSize',3);
errorbar(DfbG, minR(:,2), minCI(:,2), '--s','Color',[0.85 0.2 0.1], ...
    'LineWidth',1.4,'MarkerSize',5,'CapSize',3);
set(gca,'XScale','log');
xlabel('Position reports every $D_{\mathrm{fb}}$ intervals','Interpreter','latex');
ylabel('Minimum user rate [bit/s/Hz]','Interpreter','latex');
lg = legend({'nominal design','uncertainty-aware'},'Interpreter','latex','Location','southwest');
lg.FontSize = 7;
set(gca,'FontSize',8,'TickLabelInterpreter','latex','LineWidth',0.6, ...
    'XTick',DfbG,'XTickLabel',{'1','2','5','10','20','50'});
xlim([0.85 60]);
exportgraphics(fh, fullfile(figdir,'exp04b_sparsity.pdf'), ...
    'ContentType','vector','BackgroundColor','white','Padding','figure');

save('exp04b_feedback_sparsity.mat','DfbG','rmsePos','zPred','zSweep','minR', ...
     'minCI','pOut','outCI','itCold','itWarm','Ts','v','wNF','nEval','sig95','dPair','nEv');
fprintf('\nWrote exp04b_sparsity.pdf\n');

function [th, gb] = localMaxMin(R, th0, pOpt, N, K)
  lo = sqrt(N/K); hi = N; th = functionRobustEqualGainMM(R, th0, lo, pOpt); gb = lo;
  for it = 1:9
    g = 0.5*(lo+hi);
    t = functionRobustEqualGainMM(R, th0, g, pOpt);
    psi = zeros(K,1);
    for k = 1:K, psi(k) = sqrt(max(real(t'*R(:,:,k)*t),0)); end
    if min(psi) >= 0.98*g, lo = g; th = t; gb = g; else, hi = g; end
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
