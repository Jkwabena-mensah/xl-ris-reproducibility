addPaper2Paths();
%exp04a_rate_outage
%
%Reports the communication metrics that Section II-H of the manuscript declares
%but Section V did not carry: minimum achievable user rate, sum rate, and outage
%probability. Section V reported focusing loss and timing only, which is a gap a
%reviewer of a beam-management paper will press on.
%
%CORRECTION THIS SCRIPT CARRIES (F-35). A first version of this experiment, and
%exp03c before it, called functionRobustEqualGainMM with a FIXED target
%gbar = sqrt(N/K). That is not the algorithm of Section III-G. With a fixed
%target the equal-gain iteration is a feasibility solver: once the target is
%unreachable -- which it becomes as soon as Sigma inflates, because the moment
%kernel is sub-unit and caps every psi_k -- the objective sum_k (psi_k - gbar)^2
%degenerates into maximising sum_k psi_k, i.e. SUM gain, and the worst user is
%sacrificed. Run that way the uncertainty-aware design measures 1.4 to 3.3 dB
%WORSE than nominal, and the sum rate rises while the minimum rate falls: the
%signature of the wrong objective, not of physics. The algorithm the paper
%actually specifies wraps that iteration in a bisection on gbar, which is what
%localMaxMin below does, and under which the benefit is recovered.
%
%WHAT IS SWEPT. The dimensionless uncertainty ratio zeta = sigma_p / w_NF over
%the range for which the second-order theory of Section III is validated, with
%w_NF = lambda*rho/D formed with the aperture dimension p.D used throughout the
%paper. The traverse is zero here so that the prediction error is isolated; the
%traverse is swept in exp03g.
%
%WHAT IS COMPARED. Both designs are solved by the same max-min procedure with
%the same restarts and the same iteration budget, so the comparison isolates the
%information the design is given -- Sigma or not -- and not the solver.
%
%TWO MARGINS ARE REPORTED, DELIBERATELY.
%  criterion margin : min_k E[gain], the quantity the design optimises.
%  delivered margin : E[min_k gain], what a user actually experiences, and the
%                     quantity the rates below are built from.
%They are not the same number and the paper must not quote the first as if it
%were the second.
%
%LINK BUDGET. Absolute SNR would tie the result to a deployment. The reference is
%fixed internally: SNRop is the SNR delivered to the median user by the NOMINAL
%design at the smallest zeta on the grid. Every reported rate is then a statement
%about degradation from that operating point and is invariant to path loss,
%transmit power and noise figure, in the same spirit as the dimensionless budget M.
%
%This is version 2.0 (Last edited: 2026-09-12)
%License: GPLv2.

clear; close all; rng(11,'twister');
fprintf('exp04a v2: rate and outage, max-min by bisection\n');

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
fprintf('  N = %d, K = %d, D = %.2f m, rho = %.1f m, w_NF = %.4f m\n', N, K, p.D, rho, wNF);

zetaG = [0.02 0.04 0.06 0.08 0.10 0.12];
nZ    = numel(zetaG);
nEval = 3000;
nRest = 2;
Rth   = 2;          % outage threshold [bit/s/Hz]
SNRop = 10;         % [dB] median-user operating SNR of the reference point

gammaAll = zeros(nEval, K, nZ, 2);
critMin  = zeros(nZ,2);

tAll = tic;
for iz = 1:nZ
    sig = zetaG(iz)*wNF;  Sig = sig^2*eye(3);
    Rn = zeros(N,N,K);  Rr = zeros(N,N,K);
    for k = 1:K
        Rn(:,:,k) = functionRobustMoment(P_ris, phat(:,k), zeros(3), lam);
        Rr(:,:,k) = functionRobustMoment(P_ris, phat(:,k), Sig,      lam);
    end

    bn = -inf; br = -inf; an = []; ar = [];
    for r = 1:nRest
        th0 = exp(1j*2*pi*rand(N,1));
        tn  = localMaxMin(Rn, th0, pOpt, N, K);
        tr  = localMaxMin(Rr, th0, pOpt, N, K);
        en = inf; er = inf;
        for k = 1:K                       % scored on the TRUE moment, both designs
            en = min(en, real(tn'*Rr(:,:,k)*tn)/N^2);
            er = min(er, real(tr'*Rr(:,:,k)*tr)/N^2);
        end
        if en > bn, bn = en; an = tn; end
        if er > br, br = er; ar = tr; end
    end
    critMin(iz,:) = [bn br];

    for m = 1:nEval
        At = functionArrayResponse(P_ris, phat + sig*randn(3,K), lam);
        gammaAll(m,:,iz,1) = abs(At'*an).^2.'/N^2;
        gammaAll(m,:,iz,2) = abs(At'*ar).^2.'/N^2;
    end
    fprintf('  zeta %.2f | sigma_p %5.1f mm | criterion %+5.2f dB | %.0f s\n', ...
        zetaG(iz), 1e3*sig, 10*log10(br/bn), toc(tAll));
end

%% ---- link budget referenced to the nominal design at the smallest zeta
gRef  = median(reshape(gammaAll(:,:,1,1), [], 1));
scale = 10^(SNRop/10)/gRef;

R     = log2(1 + scale*gammaAll);
minPD = min(R, [], 2);                              % per draw, per zeta, per design
minR  = squeeze(mean(minPD, 1));
minCI = squeeze(1.96*std(minPD, 0, 1)/sqrt(nEval));
sumR  = squeeze(mean(sum(R, 2), 1));
pOut  = squeeze(mean(mean(R < Rth, 2), 1));
outCI = 1.96*sqrt(max(pOut.*(1-pOut),0)/(nEval*K));
delivdB = 10*log10(squeeze(mean(min(gammaAll,[],2),1)));

fprintf('\n zeta | criterion dB | delivered dB | min rate nom/aware   | sum rate      | P_out nom/aware\n');
for iz = 1:nZ
    fprintf(' %.2f |   %+6.2f     |   %+6.2f     | %5.2f / %5.2f +-%.3f | %5.1f / %5.1f | %.4f / %.4f\n', ...
        zetaG(iz), 10*log10(critMin(iz,2)/critMin(iz,1)), ...
        delivdB(iz,2)-delivdB(iz,1), ...
        minR(iz,1), minR(iz,2), mean(minCI(iz,:)), ...
        sumR(iz,1), sumR(iz,2), pOut(iz,1), pOut(iz,2));
end
sig95 = abs(minR(:,2)-minR(:,1)) > sqrt(minCI(:,1).^2 + minCI(:,2).^2);
fprintf('\n  min-rate difference exceeds its 95%% CI at %d of %d operating points\n', ...
    sum(sig95), nZ);

%% ---- figures
figdir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'06_Figures');
if ~isfolder(figdir), mkdir(figdir); end

fh = figure('Units','centimeters','Position',[2 2 8.8 6.8],'Color','w');
hold on; box on; grid on;
errorbar(zetaG, minR(:,1), minCI(:,1), '-o', 'Color',[0 0 0], ...
    'LineWidth',1.4,'MarkerSize',5,'CapSize',3);
errorbar(zetaG, minR(:,2), minCI(:,2), '--s', 'Color',[0.85 0.2 0.1], ...
    'LineWidth',1.4,'MarkerSize',5,'CapSize',3);
xlabel('Uncertainty ratio $\zeta = \sigma_p/w_{\mathrm{NF}}$','Interpreter','latex');
ylabel('Minimum user rate [bit/s/Hz]','Interpreter','latex');
lg = legend({'nominal design','uncertainty-aware'},'Interpreter','latex','Location','southwest');
lg.FontSize = 7;
set(gca,'FontSize',8,'TickLabelInterpreter','latex','LineWidth',0.6);
exportgraphics(fh, fullfile(figdir,'exp04a_minrate.pdf'), ...
    'ContentType','vector','BackgroundColor','white','Padding','figure');

fh2 = figure('Units','centimeters','Position',[2 2 8.8 6.8],'Color','w');
hold on; box on; grid on;
errorbar(zetaG, pOut(:,1), outCI(:,1), '-o', 'Color',[0 0 0],'LineWidth',1.4,'MarkerSize',5,'CapSize',3);
errorbar(zetaG, pOut(:,2), outCI(:,2), '--s','Color',[0.85 0.2 0.1],'LineWidth',1.4,'MarkerSize',5,'CapSize',3);
set(gca,'YScale','log');
xlabel('Uncertainty ratio $\zeta = \sigma_p/w_{\mathrm{NF}}$','Interpreter','latex');
ylabel(sprintf('Outage probability, $R_{\\mathrm{th}} = %g$ bit/s/Hz', Rth),'Interpreter','latex');
lg2 = legend({'nominal design','uncertainty-aware'},'Interpreter','latex','Location','southeast');
lg2.FontSize = 7;
set(gca,'FontSize',8,'TickLabelInterpreter','latex','LineWidth',0.6);
exportgraphics(fh2, fullfile(figdir,'exp04a_outage.pdf'), ...
    'ContentType','vector','BackgroundColor','white','Padding','figure');

save('exp04a_rate_outage.mat','zetaG','minR','minCI','sumR','pOut','outCI', ...
     'critMin','delivdB','nEval','Rth','SNRop','gRef','wNF','rho','sig95');
fprintf('\nWrote exp04a_minrate.pdf and exp04a_outage.pdf\n');

function th = localMaxMin(R, th0, pOpt, N, K)
%Max-min by bisection on the common target gbar, as specified in Section III-G.
%A fixed gbar makes the inner iteration a feasibility test, not an optimiser.
  lo = sqrt(N/K); hi = N; th = functionRobustEqualGainMM(R, th0, lo, pOpt);
  for it = 1:9
    gb = 0.5*(lo+hi);
    t  = functionRobustEqualGainMM(R, th0, gb, pOpt);
    psi = zeros(K,1);
    for k = 1:K, psi(k) = sqrt(max(real(t'*R(:,:,k)*t),0)); end
    if min(psi) >= 0.98*gb, lo = gb; th = t; else, hi = gb; end
  end
end
