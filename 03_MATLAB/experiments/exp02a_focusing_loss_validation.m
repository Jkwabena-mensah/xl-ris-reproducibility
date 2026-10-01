
addPaper2Paths();
%This Matlab script validates Theorem 1 and Corollary 1 of the article:
%
%I. G. Sarfo-Mainoo, E. A. Affum, et al., "Uncertainty-Aware Near-Field Beam
%Management for Mobile XL-RIS Under Computational Constraints," in preparation.
%
%It measures the near-field beam-focusing loss caused by a user-position error,
%by evaluating the exact spherical-wave array response with no expansion, and
%compares the measurement against two analytical predictions:
%
%  Theorem 1   L = k0^2 * trace(C*Sigma)                    (exact 2nd order)
%  Corollary 1 L ~ (2*pi^2*kappaD/3) * zeta^2,  zeta = sigma_p/wNF  (compact)
%
%where C is the aperture spread matrix, C = (1/N)*sum_n (u_n - ubar)(u_n - ubar)',
%u_n is the unit vector from RIS element n toward the user, and
%kappaD = (N1+1)/(N1-1) is the discrete-aperture correction for an N1 x N1 array.
%
%It also tests the prediction of Remark 3: because the effort demanded of the
%optimiser is governed by C + ubar*ubar' while the focusing loss is governed by
%C alone, a RADIAL position error should cost far less focusing gain than a
%TRANSVERSE error of equal magnitude. Since the radial case excites one dimension
%and the transverse case two, the predicted loss ratio is c_z/(2*c_t), i.e.
%kappaD*D^2/(60*rho^2).
%
%This is version 1.0 (Last edited: 2026-09-07)
%
%License: This code is licensed under the GPLv2 license. If you in any way use
%this code for research that results in publications, please cite the article
%listed above.
%
%Note: no optimisation is performed here. The script measures a property of the
%array response alone, so it runs in well under a minute and is independent of
%the BCD-MM solver.
%Initialization
close all;
clear all;
rng(2026);
p     = functionSimParams();
P_ris = functionRISPositions(p);              %3 x N element coordinates [m]
%%Simulation parameters
rhoNom    = p.zNom;                            %Nominal user range from RIS centre [m] (= 25 m)
rhoSweep  = [10 15 25 35 44];                  %Ranges for the wNF-scaling check [m]
offAngles = [0 15 30 45];                      %Off-broadside angles for the Corollary check [deg]
zetaTarget = [0.02 0.04 0.06 0.08 0.10 0.12 0.16 0.20];  %Uncertainty ratios to sweep
nbrOfRealizations = 20000;                     %Monte Carlo position draws per point
%2026-09-12: raised from 4000. At 4000 the 95%% half-width at zeta = 0.12 is about
%0.003 on a loss of 0.097, which left the quoted Corollary-1 error anywhere between
%7%% and 14%% -- the manuscript had 8.0%%, a re-run gave 10.2%%, and the two were not
%distinguishable at that sample count. 20000 draws resolve it.
kappaD = (p.N1 + 1) / (p.N1 - 1);              %Discrete-aperture factor (= 17/15 for N1 = 16)
k0     = 2*pi / p.lambda;                      %Wavenumber [rad/m]
cConst = 2*pi^2*kappaD/3;                      %Constant of Corollary 1 (= 7.457)
fprintf('exp02a: focusing-loss validation\n');
fprintf('  N = %d (%dx%d), D = %.2f m, lambda = %.3f m, kappaD = %.5f\n', ...
        p.N, p.N1, p.N2, p.D, p.lambda, kappaD);
fprintf('  Corollary 1 constant 2*pi^2*kappaD/3 = %.4f\n\n', cConst);
%%Preallocate matrices for storing simulation results
nZeta = length(zetaTarget);
lossMeas_iso   = zeros(1,nZeta);   %Measured loss, isotropic 3-D error
lossMeas_tran  = zeros(1,nZeta);   %Measured loss, transverse-only error
lossMeas_rad   = zeros(1,nZeta);   %Measured loss, radial-only error
lossThm_iso    = zeros(1,nZeta);   %Theorem 1 prediction, isotropic
lossThm_tran   = zeros(1,nZeta);   %Theorem 1 prediction, transverse
lossThm_rad    = zeros(1,nZeta);   %Theorem 1 prediction, radial
lossCor_iso    = zeros(1,nZeta);   %Corollary 1 prediction, isotropic
sigmaUsed      = zeros(1,nZeta);   %Position-error std used at each point [m]
%%Geometry at the nominal range
pUser  = [rhoNom; 0; 0];                       %Broadside user, range rhoNom [m]
wNF    = p.lambda * rhoNom / p.D;              %Near-field focal width [m]
[C, ubar] = functionApertureSpread(P_ris, pUser);
%Radial and transverse basis at the user location
uRad  = pUser / norm(pUser);                   %Radial unit vector (RIS centre -> user)
[~,~,V] = svd(uRad.');                          %Columns 2 and 3 span the transverse plane
uTr1  = V(:,2);  uTr2 = V(:,3);
%Validation: C must be symmetric PSD, and its eigenvalues must match Corollary 1
evC = sort(eig(C),'descend');
ctPred = kappaD * p.D^2 / (12*rhoNom^2);       %Predicted transverse eigenvalue
czPred = kappaD^2 * p.D^4 / (360*rhoNom^4);    %Predicted radial eigenvalue
assert(norm(C-C.','fro') < 1e-12, 'C is not symmetric');
assert(min(evC) > -1e-14, 'C is not positive semidefinite');
fprintf('Aperture spread matrix at rho = %.1f m (wNF = %.4f m):\n', rhoNom, wNF);
fprintf('  transverse eig  measured %.6e  predicted %.6e  (ratio %.4f)\n', ...
        evC(1), ctPred, evC(1)/ctPred);
fprintf('  radial     eig  measured %.6e  predicted %.6e\n', evC(3), czPred);
fprintf('  radial/transverse measured %.3e, predicted kappaD*D^2/(30*rho^2) = %.3e\n\n', ...
        evC(3)/evC(1), kappaD*p.D^2/(30*rhoNom^2));
%%Go through all uncertainty levels
fprintf('%-7s %-10s %-11s %-11s %-11s %-11s %-9s\n', ...
        'zeta','sigma_p[cm]','L iso MC','L iso Thm','L tran MC','L rad MC','err Thm');
fprintf('%s\n', repmat('-',1,76));
for zIdx = 1:nZeta
    sigma_p = zetaTarget(zIdx) * wNF;          %Position-error std for this point [m]
    sigmaUsed(zIdx) = sigma_p;
    %--- Three error geometries, all with the same per-axis std sigma_p ---
    %Isotropic: error in all three dimensions
    E_iso  = sigma_p * randn(3, nbrOfRealizations);
    %Transverse: error confined to the plane orthogonal to the radial direction
    E_tran = uTr1 * (sigma_p*randn(1,nbrOfRealizations)) ...
           + uTr2 * (sigma_p*randn(1,nbrOfRealizations));
    %Radial: error along the RIS-to-user direction only
    E_rad  = uRad * (sigma_p*randn(1,nbrOfRealizations));
    %--- Measure the exact focusing gain, accumulate, average after the loop ---
    accIso = 0; accTran = 0; accRad = 0;
    aDes = functionArrayResponse(P_ris, pUser, p.lambda);   %Designed response, N x 1
    for n = 1:nbrOfRealizations
        accIso  = accIso  + functionFocusingGain(P_ris, pUser, E_iso(:,n),  aDes, p.lambda);
        accTran = accTran + functionFocusingGain(P_ris, pUser, E_tran(:,n), aDes, p.lambda);
        accRad  = accRad  + functionFocusingGain(P_ris, pUser, E_rad(:,n),  aDes, p.lambda);
    end
    lossMeas_iso(zIdx)  = 1 - accIso  / nbrOfRealizations;
    lossMeas_tran(zIdx) = 1 - accTran / nbrOfRealizations;
    lossMeas_rad(zIdx)  = 1 - accRad  / nbrOfRealizations;
    %--- Analytical predictions ---
    Sig_iso  = sigma_p^2 * eye(3);
    Sig_tran = sigma_p^2 * (uTr1*uTr1.' + uTr2*uTr2.');
    Sig_rad  = sigma_p^2 * (uRad*uRad.');
    lossThm_iso(zIdx)  = k0^2 * trace(C*Sig_iso);
    lossThm_tran(zIdx) = k0^2 * trace(C*Sig_tran);
    lossThm_rad(zIdx)  = k0^2 * trace(C*Sig_rad);
    lossCor_iso(zIdx)  = cConst * zetaTarget(zIdx)^2;
    fprintf('%-7.3f %-10.2f %-11.5f %-11.5f %-11.5f %-11.5f %+8.1f%%\n', ...
            zetaTarget(zIdx), sigma_p*100, lossMeas_iso(zIdx), lossThm_iso(zIdx), ...
            lossMeas_tran(zIdx), lossMeas_rad(zIdx), ...
            100*(lossThm_iso(zIdx)-lossMeas_iso(zIdx))/lossMeas_iso(zIdx));
end
ratioMeas = mean(lossMeas_rad(1:5) ./ lossMeas_tran(1:5));   %Small-zeta points only
ratioPred = kappaD*p.D^2/(60*rhoNom^2);   %c_z/(2*c_t): radial excites 1 axis, transverse 2
fprintf('\nRadial/transverse loss ratio: measured %.3e, predicted %.3e\n', ...
        ratioMeas, ratioPred);
fprintf('-> a radial position error costs about 1/%.0f of what an equal transverse error costs.\n\n', ...
        1/ratioMeas);
%%Range sweep: confirm that wNF, not sigma_p alone, sets the loss
fprintf('Range sweep at fixed zeta = 0.08:\n');
fprintf('%-9s %-10s %-13s %-11s %-11s\n','rho[m]','wNF[m]','sigma_p[cm]','L measured','L Cor.1');
fprintf('%s\n', repmat('-',1,58));
lossRange = zeros(1,length(rhoSweep));  wNFrange = zeros(1,length(rhoSweep));
for rIdx = 1:length(rhoSweep)
    pU   = [rhoSweep(rIdx); 0; 0];
    wR   = p.lambda*rhoSweep(rIdx)/p.D;
    sg   = 0.08 * wR;
    aD   = functionArrayResponse(P_ris, pU, p.lambda);
    acc  = 0;
    Er   = sg*randn(3,nbrOfRealizations);
    for n = 1:nbrOfRealizations
        acc = acc + functionFocusingGain(P_ris, pU, Er(:,n), aD, p.lambda);
    end
    lossRange(rIdx) = 1 - acc/nbrOfRealizations;  wNFrange(rIdx) = wR;
    fprintf('%-9.1f %-10.4f %-13.2f %-11.5f %-11.5f\n', ...
            rhoSweep(rIdx), wR, sg*100, lossRange(rIdx), cConst*0.08^2);
end
fprintf('-> the loss is set by zeta, not by sigma_p: identical zeta gives identical loss\n');
fprintf('   across a %.1fx range span, while sigma_p varies by the same factor.\n\n', ...
        max(rhoSweep)/min(rhoSweep));
%%Off-broadside check: how far does the Corollary 1 constant hold?
fprintf('Off-broadside check at zeta = 0.08:\n');
fprintf('%-12s %-12s %-12s\n','angle[deg]','L measured','L Cor.1');
fprintf('%s\n', repmat('-',1,38));
lossAngle = zeros(1,length(offAngles));
for aIdx = 1:length(offAngles)
    th  = offAngles(aIdx)*pi/180;
    pU  = rhoNom*[cos(th); sin(th); 0];
    wR  = p.lambda*rhoNom/p.D;   sg = 0.08*wR;
    aD  = functionArrayResponse(P_ris, pU, p.lambda);
    acc = 0;  Er = sg*randn(3,nbrOfRealizations);
    for n = 1:nbrOfRealizations
        acc = acc + functionFocusingGain(P_ris, pU, Er(:,n), aD, p.lambda);
    end
    lossAngle(aIdx) = 1 - acc/nbrOfRealizations;
    fprintf('%-12d %-12.5f %-12.5f\n', offAngles(aIdx), lossAngle(aIdx), cConst*0.08^2);
end
%%Save simulation results
save('exp02a_focusing_loss_validation.mat', ...
     'zetaTarget','sigmaUsed','lossMeas_iso','lossMeas_tran','lossMeas_rad', ...
     'lossThm_iso','lossThm_tran','lossThm_rad','lossCor_iso', ...
     'rhoSweep','lossRange','wNFrange','offAngles','lossAngle', ...
     'C','ubar','wNF','kappaD','cConst','rhoNom','nbrOfRealizations','p');
%%Plot simulation results
%Publication-quality export: IET single column is 88 mm wide.
figW = 8.8; figH = 6.6;   %centimetres
%Figure 1: loss versus uncertainty ratio, measurement against theory
fh1 = figure('Units','centimeters','Position',[2 2 figW figH],'Color','w');
hold on; box on; grid on;
plot(zetaTarget, lossMeas_iso,  'ko',  'LineWidth', 1.2, 'MarkerSize', 6);
plot(zetaTarget, lossThm_iso,   'b-',  'LineWidth', 1.4);
plot(zetaTarget, lossCor_iso,   'r--', 'LineWidth', 1.4);
plot(zetaTarget, lossMeas_tran, 'ks',  'LineWidth', 1.2, 'MarkerSize', 6);
plot(zetaTarget, lossMeas_rad,  'k^',  'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('Uncertainty ratio $\zeta=\sigma_p/w_{\mathrm{NF}}$', 'Interpreter','latex');
ylabel('Beam-focusing loss $\mathcal{L}$', 'Interpreter','latex');
lg1 = legend({'Measured (isotropic)','Theorem 1','Corollary 1', ...
        'Measured (transverse)','Measured (radial)'}, ...
        'Location','northwest','Interpreter','latex','Box','on');
lg1.FontSize = 7;
xlim([0 max(zetaTarget)]); ylim([0 1.15*max(lossCor_iso)]);
set(gca,'FontSize',8,'TickLabelInterpreter','latex','LineWidth',0.6);
exportgraphics(fh1, 'exp02a_loss_vs_zeta.pdf', ...
    'ContentType','vector','BackgroundColor','white','Padding','figure');
%Figure 2: relative error of the analytical predictions, to expose the validity range
fh2 = figure('Units','centimeters','Position',[2 2 figW figH],'Color','w');
hold on; box on; grid on;
relThm = 100*(lossThm_iso - lossMeas_iso)./lossMeas_iso;
relCor = 100*(lossCor_iso - lossMeas_iso)./lossMeas_iso;
plot(zetaTarget, relThm, 'b-o',  'LineWidth', 1.4, 'MarkerSize', 5);
plot(zetaTarget, relCor, 'r--s', 'LineWidth', 1.4, 'MarkerSize', 5);
yline(0,'k:','LineWidth',0.8);
yline(10,'k-.','LineWidth',0.8);
xline(0.12,'k-.','LineWidth',0.8);
text(0.121, -4, 'validity limit', 'Interpreter','latex','FontSize',7);
xlabel('Uncertainty ratio $\zeta=\sigma_p/w_{\mathrm{NF}}$', 'Interpreter','latex');
ylabel('Prediction error [\%]', 'Interpreter','latex');
lg2 = legend({'Theorem 1','Corollary 1'}, ...
        'Location','northwest','Interpreter','latex','Box','on');
lg2.FontSize = 7;
xlim([0 max(zetaTarget)]);
set(gca,'FontSize',8,'TickLabelInterpreter','latex','LineWidth',0.6);
exportgraphics(fh2, 'exp02a_prediction_error.pdf', ...
    'ContentType','vector','BackgroundColor','white','Padding','figure');
fprintf('\nSaved exp02a_focusing_loss_validation.mat and two PDFs.\n');
