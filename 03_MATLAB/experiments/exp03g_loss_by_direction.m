%exp03g_loss_by_direction
%
%Dedicated, high-sample computation of the focusing loss over the interval as a
%function of reconfiguration interval and direction of travel. Split out from
%exp03f because that experiment evaluated the loss with 120 draws per slot, which
%left enough Monte-Carlo noise for the transverse curve to run non-monotonically
%in Ts -- an artefact, not physics, and not something to put in a figure. Here the
%designs are unchanged but the evaluation uses 3000 draws and reports a 95%
%confidence interval, so that what is plotted can be distinguished from noise.
%
%This is version 1.0 (Last edited: 2026-09-11)
%License: GPLv2.

clear; close all; rng(31,'twister');
addPaper2Paths();
fprintf('exp03g: focusing loss by direction, high sample count\n');

p = functionSimParams(); P_ris = functionRISPositions(p);
bin = functionArrayResponse(P_ris, p.pBS, p.lambda);
N = p.N; K = p.K; lam = p.lambda; gbar = sqrt(N/K);
pOpt = p; pOpt.Imax = 400; pOpt.muTol = 1e-4;

pU0 = zeros(3,K);
for k = 1:K
    pU0(:,k) = [p.region(1,1)+diff(p.region(1,:))*rand;
                p.region(2,1)+diff(p.region(2,:))*rand;
                p.region(3,1)+diff(p.region(3,:))*rand];
end
th0 = exp(1j*2*pi*rand(N,1));

TsG  = [0.5 1 2 3 5]*1e-3;
dirs = {[1;0;0],'radial'; [1;1;0]/sqrt(2),'45deg'; [0;1;0],'transverse'};
v = 60/3.6; Dfb = 1; nSlot = 5; nEval = 3000;
L = zeros(numel(TsG),3); Lci = zeros(numel(TsG),3);

for it = 1:numel(TsG)
  Ts = TsG(it);
  Sig = localCov(Ts,Dfb,1.0,p.r); sigp = sqrt(trace(Sig)/3);
  for id = 1:3
    vh = dirs{id,1}; delta = v*Ts*Dfb;
    pk = pU0; th = th0; acc = [];
    for s = 1:nSlot
      pk = pk + vh*(v*Ts);
      Fk = cell(1,K);
      for k = 1:K
        Fk{k} = functionRobustMomentFactor(P_ris,pk(:,k),zeros(3),lam,0);
      end
      th = functionRobustEqualGainMM(Fk,th,gbar,pOpt);
      if s > 2
        gk = zeros(nEval,K);
        for m = 1:nEval
          e = sigp*randn(3,K) + vh*((rand(1,K)-0.5)*delta);
          At = functionArrayResponse(P_ris, pk+e, lam);
          gk(m,:) = abs(At'*th).^2.'/N^2;
        end
        g0 = abs(functionArrayResponse(P_ris,pk,lam)'*th).^2.'/N^2;
        lk = max(0, 1 - mean(gk,1)./max(g0,eps));      %per-user loss
        [lw, kw] = max(lk);                            %worst user
        se = std(gk(:,kw))/sqrt(nEval)/max(g0(kw),eps);
        acc(end+1,:) = [lw, 1.96*se];                  %#ok<AGROW>
      end
    end
    L(it,id)   = mean(acc(:,1));
    Lci(it,id) = mean(acc(:,2));
  end
  fprintf('  Ts = %4.1f ms | radial %5.2f%%  45deg %5.2f%%  transverse %5.2f%%  (+-%.2f)\n', ...
      Ts*1e3, 100*L(it,1), 100*L(it,2), 100*L(it,3), 100*mean(Lci(it,:)));
end

%% figure
figdir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'06_Figures');
if ~isfolder(figdir), mkdir(figdir); end
fh = figure('Units','centimeters','Position',[2 2 8.8 6.8],'Color','w');
hold on; box on; grid on;
mk = {'-o','-s','-^'}; col = [0 0 0; 0 0.35 0.85; 0.85 0.2 0.1];
for id = 1:3
    errorbar(TsG*1e3, 100*L(:,id), 100*Lci(:,id), mk{id}, ...
        'Color', col(id,:), 'LineWidth',1.4, 'MarkerSize',5, 'CapSize',3);
end
yline(10,'k--','LineWidth',1.0);
text(3.05, 11.4, 'loss budget','Interpreter','latex','FontSize',7);
set(gca,'YScale','log');
xlabel('Reconfiguration interval $T_s$ [ms]','Interpreter','latex');
ylabel('Focusing loss over the interval [\%]','Interpreter','latex');
lg = legend({'radial','$45^\circ$','transverse'},'Interpreter','latex','Location','northwest');
lg.FontSize = 7;
xlim([0.4 5.4]); ylim([0.2 20]);
set(gca,'FontSize',8,'TickLabelInterpreter','latex','LineWidth',0.6, ...
    'XTick',[0.5 1 2 3 5],'YTick',[0.3 1 3 10]);
exportgraphics(fh, fullfile(figdir,'exp03g_loss_direction.pdf'), ...
    'ContentType','vector','BackgroundColor','white','Padding','figure');
save('exp03g_loss.mat','TsG','dirs','L','Lci','nEval','v');
fprintf('\nWrote exp03g_loss_direction.pdf\n');

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
