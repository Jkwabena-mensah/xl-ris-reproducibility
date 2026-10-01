%plot_ieee_figs_v9  Regenerate the whole Paper 2 figure set in one pass.
%
%WHAT THIS DOES AND DOES NOT DO. Every figure below is drawn from the .mat file
%that the corresponding experiment already wrote. Nothing is re-simulated and no
%stored quantity is recomputed, so no number in the manuscript can move as a
%result of running this. What changes is the plot TYPE where the previous type
%could not show the result, and the visual vocabulary, which is now shared
%across all six figures.
%
%THE ENCODING RULE. No series is distinguished by colour alone. Each carries a
%colour, a line style and a marker together, so the set survives a greyscale
%print and a colour-blind reader. See functionIEEEPalette.
%
%DRAWING ORDER. Error bars are drawn first and markers on top of them. A bar
%shorter than the marker would otherwise be painted inside the marker's face and
%fill what should read as a hollow symbol -- which is what happened at
%D_fb = 50, where the confidence interval is a tenth of the marker's height.
%
%FIGURES 6, 7 AND 8 NOW COME FROM THE RESTART-SAFE EXPERIMENTS exp12a AND
%exp12b, not from exp04a and exp04b. Section V-E's prose was written from
%exp04c while its figures were generated from exp04a -- nine numbers matched the
%former and none the latter -- and both source experiments selected the best of
%R restarts under the swept-parameter moment, which is the confound that
%manufactured the Section V-G reversal. exp12a and exp12b solve the nominal bank
%once and select nothing, so text and figures finally describe one experiment.
%
%MARKER CONVENTION on every paired panel: a FILLED marker means the 95%%
%confidence interval on the paired difference excludes zero, a hollow marker
%means it does not. The verdict is in the figure, not only in the caption.
%
%THE THREE TYPE CHANGES.
%  Fig. 3  linear -> logarithmic ordinate. The radial series spans 1.3e-7 to
%          1.3e-5 while the isotropic series reaches 0.3; on a linear axis the
%          radial measurements lie on the frame and the figure cannot show the
%          direction dependence its caption claims. Four decades of separation
%          are the result.
%  Fig. 7  the two smallest uncertainty levels recorded ZERO outage events. The
%          previous figure plotted the rule-of-three bound at those points as
%          though it were a measured probability, joined by the same line. They
%          are now drawn as censored observations -- hollow marker at the bound
%          with a downward arrow -- and the line is broken there.
%  Fig. 8  gains a lower panel carrying the paired difference that exp04b
%          already computes and stores in dPair. The upper panel's error bars
%          are the spread of each design's absolute rate, which is far wider
%          than the paired contrast and hides a separation that is in fact
%          decisive at large feedback intervals.
%
%WHAT CHANGED IN v9. Figures 6, 7, 8 and 9 now come from exp14a, exp14b and
%exp14e, which differ from exp12a/exp12b/exp11b in two respects: the max-min
%bisection returns the best probe rather than the last accepted one, and the
%restart bank is 24 rather than 12, and 48 for Figure 9. Figure 7 gains the worst-case ellipsoidal
%design of Section V-H as a third series, and its censoring is now handled per
%series because the three designs fall silent at different uncertainty levels.
%
%v11 applies the IEEE figure checklist: units in parentheses rather than square
%brackets, legend type raised to the 8 pt floor with the widest labels shortened
%to fit, Fig. 3's legend moved inside the axes (northoutside was spending 17%% of
%a 7.2 cm figure on a legend while the axes held a large empty band), and Fig. 7's
%ordinate reduced to the quantity plotted, its threshold belonging to the caption.
%No data, no series, no limits and no export size change.
%
%The two-panel ordinates are set on two lines: parenthesised units are wider than
%the bracketed form they replace, and the single-line label overran its panel.
%
%Version 1.5 (2026-09-19). License: GPLv2.

clear; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % package root (was a hard-coded C:\Dev path; release edit)
expd = fullfile(root,'03_MATLAB','experiments');
figd = fullfile(root,'07_Manuscript','figures');
addpath(expd);
cd(expd);
S = functionIEEEPalette();
fprintf('\n=== regenerating the Paper 2 figure set ===\n');

%% ===== Fig. 3 -- focusing loss against the uncertainty ratio ==============
A = load('exp02a_focusing_loss_validation.mat');
z = A.zetaTarget(:);

fh = figure('Units','centimeters','Position',[2 2 S.wCol 7.2],'Color','w');
ax = axes(fh); hold(ax,'on');
set(ax,'YScale','log');

%Theorem 1 predicts the transverse case to within 4.3e-5 relative of the
%isotropic case, so drawing both lines would lay one exactly over the other.
%The prediction is drawn once and both sets of measurements are shown against it.
plot(ax,z,A.lossThm_iso,S.s1.ls,'Color',S.s1.c,'LineWidth',S.lw,'HandleVisibility','off');
plot(ax,z,A.lossThm_rad,S.s3.ls,'Color',S.s3.c,'LineWidth',S.lw,'HandleVisibility','off');

%Corollary 1 over the Theorem 1 line it reproduces, so its dots read through
plot(ax,z,A.lossCor_iso,S.s4.ls,'Color',S.s4.c,'LineWidth',1.4,'HandleVisibility','off');

%simulation (markers only, drawn last so nothing paints over their faces)
plot(ax,z,A.lossMeas_iso ,S.s1.mk,'Color',S.s1.c,'MarkerFaceColor','w', ...
     'MarkerSize',S.ms,'LineWidth',0.9,'LineStyle','none','HandleVisibility','off');
plot(ax,z,A.lossMeas_tran,S.s2.mk,'Color',S.s2.c,'MarkerFaceColor','w', ...
     'MarkerSize',S.ms,'LineWidth',0.9,'LineStyle','none','HandleVisibility','off');
plot(ax,z,A.lossMeas_rad ,S.s3.mk,'Color',S.s3.c,'MarkerFaceColor','w', ...
     'MarkerSize',S.ms,'LineWidth',0.9,'LineStyle','none','HandleVisibility','off');

%legend proxies: line = Theorem 1, marker = Monte Carlo, colour = direction
p1 = plot(ax,nan,nan,[S.s1.ls S.s1.mk],'Color',S.s1.c,'MarkerFaceColor','w', ...
          'LineWidth',S.lw,'MarkerSize',S.ms);
p2 = plot(ax,nan,nan,S.s2.mk,'Color',S.s2.c,'MarkerFaceColor','w', ...
          'LineStyle','none','MarkerSize',S.ms);
p3 = plot(ax,nan,nan,[S.s3.ls S.s3.mk],'Color',S.s3.c,'MarkerFaceColor','w', ...
          'LineWidth',S.lw,'MarkerSize',S.ms);
p4 = plot(ax,nan,nan,S.s4.ls,'Color',S.s4.c,'LineWidth',1.4);

xlabel(ax,'Uncertainty ratio $\zeta=\sigma_p/w_{\mathrm{NF}}$','Interpreter','latex');
ylabel(ax,'Beam-focusing loss $\mathcal{L}$','Interpreter','latex');
xlim(ax,[0 0.21]); xticks(ax,0:0.05:0.20); xticklabels(ax,compose('%.2f',(0:0.05:0.20)'));
ylim(ax,[3e-8 1]); yticks(ax,10.^(-7:0));
legend(ax,[p1 p2 p3 p4],{'Isotropic','Transverse','Radial','Corollary 1'}, ...
       'Interpreter','latex','Location','east','NumColumns',1);
functionIEEEFinish(fh, fullfile(figd,'exp02a_loss_vs_zeta.pdf'), S.wCol, 7.2);

%% ===== Fig. 4 -- accuracy of the closed-form prediction ==================
relThm = 100*(A.lossThm_iso./A.lossMeas_iso - 1);
relCor = 100*(A.lossCor_iso./A.lossMeas_iso - 1);

fh = figure('Units','centimeters','Position',[2 2 S.wCol 6.2],'Color','w');
ax = axes(fh); hold(ax,'on');

%the validated operating range, shaded rather than marked by a bare rule
fill(ax,[0 0.12 0.12 0],[-6 -6 32 32],S.teal,'FaceAlpha',0.075, ...
     'EdgeColor','none','HandleVisibility','off');
yline(ax,0 ,'-' ,'Color',S.grey,'LineWidth',0.6,'HandleVisibility','off');
yline(ax,10,'--','Color',S.grey,'LineWidth',0.8,'HandleVisibility','off');

h1 = localSeries(ax, z, relThm, [], [], S.s1, S);
h2 = localSeries(ax, z, relCor, [], [], S.s4, S);

text(ax,0.060,-3.4,'validated range, $\zeta\le0.12$','Interpreter','latex', ...
     'FontSize',S.fsAnn,'Color',S.teal*0.8,'HorizontalAlignment','center');
text(ax,0.205,10,'$10\%$','Interpreter','latex','FontSize',S.fsAnn, ...
     'Color',S.grey*0.7,'HorizontalAlignment','right','VerticalAlignment','bottom');

xlabel(ax,'Uncertainty ratio $\zeta=\sigma_p/w_{\mathrm{NF}}$','Interpreter','latex');
ylabel(ax,'Prediction error (\%)','Interpreter','latex');
xlim(ax,[0 0.21]); xticks(ax,0:0.05:0.20); xticklabels(ax,compose('%.2f',(0:0.05:0.20)'));
ylim(ax,[-6 32]); yticks(ax,-5:5:30);
legend(ax,[h1 h2],{'Theorem 1','Corollary 1'},'Interpreter','latex', ...
       'Location','northwest');
functionIEEEFinish(fh, fullfile(figd,'exp02a_prediction_error.pdf'), S.wCol, 6.2);

%% ===== Fig. 6 -- minimum user rate, absolute and paired =================
A  = load('exp14a_bestiter.mat');
zb = A.zetaG(:);
mA = mean(A.dAN,2); sA = 1.96*std(A.dAN,0,2)/sqrt(A.RB);
okA = (mA-sA > 0) | (mA+sA < 0);              % paired interval excludes zero
rN = mean(A.rNom,2); rA = mean(A.rAwa,2);
cN = 1.96*std(A.rNom,0,2)/sqrt(A.RB); cA = 1.96*std(A.rAwa,0,2)/sqrt(A.RB);

fh = figure('Units','centimeters','Position',[2 2 S.wCol 8.2],'Color','w');
tl = tiledlayout(fh,2,1,'TileSpacing','compact','Padding','compact');

axa = nexttile(tl); hold(axa,'on');
u1 = localSeries(axa, zb, rN, cN, cN, S.s1, S);
u2 = localSeries(axa, zb, rA, cA, cA, S.s2, S);
ylabel(axa,{'Minimum user rate','(bit/s/Hz)'},'Interpreter','latex');
xlim(axa,[0.012 0.128]); xticks(axa,zb); xticklabels(axa,{});
ylim(axa,[2.58 3.66]); yticks(axa,2.75:0.25:3.5);
legend(axa,[u1 u2],{'Nominal design','Uncertainty-aware'},'Interpreter','latex', ...
       'Location','southwest');
text(axa,0.930,0.86,'(a)','Units','normalized','Interpreter','latex', ...
     'FontSize',S.fsLab,'FontName',S.fnt,'HorizontalAlignment','right');

axb = nexttile(tl); hold(axb,'on');
yline(axb,0,'-','Color',S.grey,'LineWidth',0.6,'HandleVisibility','off');
localPaired(axb, zb, mA, sA, okA, S.s3, S);
xlabel(axb,'Uncertainty ratio $\zeta=\sigma_p/w_{\mathrm{NF}}$','Interpreter','latex');
ylabel(axb,{'Paired gain in','minimum rate (bit/s/Hz)'},'Interpreter','latex');
xlim(axb,[0.012 0.128]); xticks(axb,zb); xticklabels(axb,compose('%.2f',zb));
ylim(axb,[-0.012 0.096]); yticks(axb,0:0.02:0.08);
text(axb,0.930,0.86,'(b)','Units','normalized','Interpreter','latex', ...
     'FontSize',S.fsLab,'FontName',S.fnt,'HorizontalAlignment','right');
functionIEEEFinish(fh, fullfile(figd,'exp04a_minrate.pdf'), S.wCol, 8.2);

%% ===== Fig. 7 -- outage, three designs, censoring handled per series =====
%THE CENSORING IS PER SERIES. The nominal design records its first event at
%zeta = 0.04, the uncertainty-aware design at 0.06 and the worst-case design at
%0.08. Plotting a common floor would assert three measurements where there are
%none; each series is therefore broken where it recorded nothing, and the
%silence is drawn as a hollow marker at the rule-of-three bound with a downward
%arrow.
nPool  = A.nEval*A.K*A.RB;             % 24 restarts x 24 000 user-slots
floor3 = 3/nPool;                      % rule-of-three 95% bound for zero events
P3     = [mean(A.oNom,2), mean(A.oAwa,2), mean(A.oWc,2)];
spec   = {S.s1, S.s2, S.s4};
dxs    = [-1 0 1]*0.0016;

fh = figure('Units','centimeters','Position',[2 2 S.wCol 6.2],'Color','w');
ax = axes(fh); hold(ax,'on');
set(ax,'YScale','log');
yline(ax,floor3,':','Color',S.grey,'LineWidth',0.9,'HandleVisibility','off');

hh = gobjects(1,3);
for d = 1:3
    msr = P3(:,d) > 0;
    hh(d) = localSeries(ax, zb(msr), P3(msr,d), [], [], spec{d}, S);
    for ii = find(~msr)'
        xx = zb(ii) + dxs(d);
        plot(ax,[xx xx],[floor3 floor3/3.2],'-','Color',spec{d}.c, ...
             'LineWidth',0.8,'HandleVisibility','off');
        plot(ax,xx,floor3/3.2,'v','Color',spec{d}.c,'MarkerFaceColor',spec{d}.c, ...
             'MarkerSize',3,'LineStyle','none','HandleVisibility','off');
        plot(ax,xx,floor3,spec{d}.mk,'Color',spec{d}.c,'MarkerFaceColor','w', ...
             'MarkerSize',S.ms,'LineWidth',0.9,'LineStyle','none', ...
             'HandleVisibility','off');
    end
end
h4 = plot(ax,nan,nan,'v','Color',S.grey,'MarkerFaceColor',S.grey,'MarkerSize',3, ...
          'LineStyle','none');
text(ax,0.126,floor3/1.45,'rule-of-three bound, $n=576\,000$', ...
     'Interpreter','latex','FontSize',S.fsAnn,'Color',S.grey*0.7, ...
     'HorizontalAlignment','right','VerticalAlignment','top');
xlabel(ax,'Uncertainty ratio $\zeta=\sigma_p/w_{\mathrm{NF}}$','Interpreter','latex');
ylabel(ax,'Outage probability','Interpreter','latex');
xlim(ax,[0.012 0.128]); xticks(ax,zb); xticklabels(ax,compose('%.2f',zb));
ylim(ax,[6e-7 6e-2]); yticks(ax,10.^(-6:-2));
legend(ax,[hh h4],{'Nominal design','Uncertainty-aware','Worst case, $\mathcal{U}_k$', ...
       'No event observed'},'Interpreter','latex','Location','northwest');
functionIEEEFinish(fh, fullfile(figd,'exp04a_outage.pdf'), S.wCol, 6.2);

%% ===== Fig. 8 -- feedback sparsity, absolute and paired =================
C  = load('exp14b_bestiter.mat');
D  = C.DfbG(:);
mC = mean(C.dRate,2); sC = 1.96*std(C.dRate,0,2)/sqrt(C.RB);
okC = (mC-sC > 0) | (mC+sC < 0);
qN = mean(C.rNom,2); qA = mean(C.rAwa,2);
eN = 1.96*std(C.rNom,0,2)/sqrt(C.RB); eA = 1.96*std(C.rAwa,0,2)/sqrt(C.RB);

fh = figure('Units','centimeters','Position',[2 2 S.wCol 8.2],'Color','w');
tl = tiledlayout(fh,2,1,'TileSpacing','compact','Padding','compact');

axa = nexttile(tl); hold(axa,'on'); set(axa,'XScale','log');
q1 = localSeries(axa, D, qN, eN, eN, S.s1, S);
q2 = localSeries(axa, D, qA, eA, eA, S.s2, S);
ylabel(axa,{'Minimum user rate','(bit/s/Hz)'},'Interpreter','latex');
xlim(axa,[0.85 60]); xticks(axa,D); xticklabels(axa,{});
ylim(axa,[0.5 4.05]); yticks(axa,1:0.5:3.5);
legend(axa,[q1 q2],{'Nominal design','Uncertainty-aware'},'Interpreter','latex', ...
       'Location','southwest');
text(axa,0.030,0.90,'(a)','Units','normalized','Interpreter','latex', ...
     'FontSize',S.fsLab,'FontName',S.fnt);

axb = nexttile(tl); hold(axb,'on'); set(axb,'XScale','log');
yline(axb,0,'-','Color',S.grey,'LineWidth',0.6,'HandleVisibility','off');
localPaired(axb, D, mC, sC, okC, S.s3, S);
xlabel(axb,'Position reports every $D_{\mathrm{fb}}$ intervals','Interpreter','latex');
ylabel(axb,{'Paired gain in','minimum rate (bit/s/Hz)'},'Interpreter','latex');
xlim(axb,[0.85 60]); xticks(axb,D); xticklabels(axb,compose('%d',D));
ylim(axb,[-0.08 0.58]); yticks(axb,0:0.1:0.5);
text(axb,0.030,0.90,'(b)','Units','normalized','Interpreter','latex', ...
     'FontSize',S.fsLab,'FontName',S.fnt);
functionIEEEFinish(fh, fullfile(figd,'exp04b_sparsity.pdf'), S.wCol, 8.2);

%% ===== Fig. 9 -- is the comparison resolved, and under which model? ====
%PLOT TYPE CHANGED. The claim this figure carries is not about two means; it is
%about which of two comparisons is RESOLVED against solver restart variability.
%A line-and-whisker plot invites the reader to compare the means and leaves the
%verdict to the caption. Panel (a) therefore shows every restart as a point, so
%the spread is visible as data rather than as a whisker, and panel (b) plots the
%statistic the verdict actually rests on against the band in which nothing is
%resolved. The line-of-sight series leaves that band at zeta = 0.04 and stays
%out; the cascaded-Rician series never leaves it.
E  = load('exp14e_v5g48.mat');
ze = E.zetaG(:);
mL = mean(E.VL,2); sL = 1.96*std(E.VL,0,2)/sqrt(E.RB); tL = mL./(sL/1.96);
mR = mean(E.VR,2); sR = 1.96*std(E.VR,0,2)/sqrt(E.RB); tR = mR./(sR/1.96);
de = 0.0022;                        % dodge between the two evaluation models
rng(9,'twister'); jit = @(n) (rand(n,1)-0.5)*0.0018;   % fixed, so the figure redraws identically

fh = figure('Units','centimeters','Position',[2 2 S.wCol 8.2],'Color','w');
tl = tiledlayout(fh,2,1,'TileSpacing','compact','Padding','compact');

%---- (a) every restart, with the mean and its 95% interval on top
axa = nexttile(tl); hold(axa,'on');
yline(axa,0,'-','Color',S.grey,'LineWidth',0.6,'HandleVisibility','off');
for i = 1:numel(ze)
    scatter(axa, ze(i)-de+jit(E.RB), E.VL(i,:)', 8, 'o', ...
        'MarkerEdgeColor','none','MarkerFaceColor',S.s1.c,'MarkerFaceAlpha',0.38, ...
        'HandleVisibility','off');
    scatter(axa, ze(i)+de+jit(E.RB), E.VR(i,:)', 8, 's', ...
        'MarkerEdgeColor','none','MarkerFaceColor',S.s2.c,'MarkerFaceAlpha',0.38, ...
        'HandleVisibility','off');
end
errorbar(axa,ze+de,mR,sR,'LineStyle',S.s2.ls,'Color',S.s2.c,'Marker','none', ...
     'LineWidth',S.lw,'CapSize',S.cap,'HandleVisibility','off');
errorbar(axa,ze-de,mL,sL,'LineStyle',S.s1.ls,'Color',S.s1.c,'Marker','none', ...
     'LineWidth',S.lw,'CapSize',S.cap,'HandleVisibility','off');
plot(axa,ze+de,mR,S.s2.mk,'Color',S.s2.c,'MarkerFaceColor','w','MarkerSize',S.ms, ...
     'LineWidth',0.9,'LineStyle','none','HandleVisibility','off');
plot(axa,ze-de,mL,S.s1.mk,'Color',S.s1.c,'MarkerFaceColor','w','MarkerSize',S.ms, ...
     'LineWidth',0.9,'LineStyle','none','HandleVisibility','off');
g1 = plot(axa,nan,nan,[S.s1.ls S.s1.mk],'Color',S.s1.c,'MarkerFaceColor','w', ...
     'LineWidth',S.lw,'MarkerSize',S.ms);
g2 = plot(axa,nan,nan,[S.s2.ls S.s2.mk],'Color',S.s2.c,'MarkerFaceColor','w', ...
     'LineWidth',S.lw,'MarkerSize',S.ms);
ylabel(axa,{'Paired gain in','minimum rate (bit/s/Hz)'},'Interpreter','latex');
xlim(axa,[0.012 0.128]); xticks(axa,ze); xticklabels(axa,{});
ylim(axa,[-1.15 1.15]); yticks(axa,-1:0.5:1);
legend(axa,[g1 g2],{'Line-of-sight','Cascaded Rician'}, ...
       'Interpreter','latex','Location','northwest','NumColumns',1);
text(axa,0.965,0.92,'(a)','Units','normalized','Interpreter','latex', ...
     'FontSize',S.fsLab,'FontName',S.fnt,'HorizontalAlignment','right');

%---- (b) the verdict: the statistic against the band where nothing is resolved
axb = nexttile(tl); hold(axb,'on');
fill(axb,[0.012 0.128 0.128 0.012],[-1.96 -1.96 1.96 1.96],S.grey, ...
     'FaceAlpha',0.16,'EdgeColor','none','HandleVisibility','off');
yline(axb,0,'-','Color',S.grey,'LineWidth',0.6,'HandleVisibility','off');
localSeries(axb, ze, tR, [], [], S.s2, S);
localSeries(axb, ze, tL, [], [], S.s1, S);
text(axb,0.070,-1.20,'not resolved ($|t|<1.96$)','Interpreter','latex', ...
     'FontSize',S.fsAnn,'Color',[0.35 0.35 0.35],'HorizontalAlignment','center');
xlabel(axb,'Uncertainty ratio $\zeta=\sigma_p/w_{\mathrm{NF}}$','Interpreter','latex');
ylabel(axb,{'Paired $t$ statistic','over 48 restarts'},'Interpreter','latex');
xlim(axb,[0.012 0.128]); xticks(axb,ze); xticklabels(axb,compose('%.2f',ze));
ylim(axb,[-1.55 6.55]); yticks(axb,-1:1:6);
text(axb,0.035,0.92,'(b)','Units','normalized','Interpreter','latex', ...
     'FontSize',S.fsLab,'FontName',S.fnt,'HorizontalAlignment','left');
functionIEEEFinish(fh, fullfile(figd,'exp04c_rician.pdf'), S.wCol, 8.2);

fprintf('\n=== done. six figures written to %s ===\n', figd);
fprintf('LoS t: %s\n', mat2str(round((mL./(sL/1.96))',2)));
fprintf('Ric t: %s\n', mat2str(round((mR./(sR/1.96))',2)));

%% ================= local functions =================
function h = localSeries(ax, x, y, eLo, eHi, spec, S)
%Draw one series as connecting line, then error bars, then markers, in that
%order, and return a legend proxy carrying the combined line-and-marker glyph.
%The order is the point: an error bar shorter than the marker would otherwise
%be painted across the marker's face and destroy the hollow symbol.
    if isempty(eLo)
        plot(ax,x,y,spec.ls,'Color',spec.c,'LineWidth',S.lw,'HandleVisibility','off');
    else
        errorbar(ax,x,y,eLo,eHi,'LineStyle',spec.ls,'Color',spec.c, ...
                 'Marker','none','LineWidth',S.lw,'CapSize',S.cap, ...
                 'HandleVisibility','off');
    end
    plot(ax,x,y,spec.mk,'Color',spec.c,'MarkerFaceColor','w','MarkerSize',S.ms, ...
         'LineWidth',0.9,'LineStyle','none','HandleVisibility','off');
    h = plot(ax,nan,nan,[spec.ls spec.mk],'Color',spec.c,'MarkerFaceColor','w', ...
             'LineWidth',S.lw,'MarkerSize',S.ms);
end

function localPaired(ax, x, m, s, ok, spec, S)
%A paired-difference series. The marker is FILLED where the 95% interval
%excludes zero and hollow where it does not, so the resolution verdict is
%carried by the figure rather than by the caption.
    errorbar(ax,x,m,s,s,'LineStyle',spec.ls,'Color',spec.c,'Marker','none', ...
             'LineWidth',S.lw,'CapSize',S.cap,'HandleVisibility','off');
    if any(~ok)
        plot(ax,x(~ok),m(~ok),spec.mk,'Color',spec.c,'MarkerFaceColor','w', ...
             'MarkerSize',S.ms,'LineWidth',0.9,'LineStyle','none','HandleVisibility','off');
    end
    if any(ok)
        plot(ax,x(ok),m(ok),spec.mk,'Color',spec.c,'MarkerFaceColor',spec.c, ...
             'MarkerSize',S.ms,'LineWidth',0.9,'LineStyle','none','HandleVisibility','off');
    end
end
