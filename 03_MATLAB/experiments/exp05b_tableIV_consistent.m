%exp05b_tableIV_consistent
%
%Table IV and Section III-F disagreed because they were measured in different
%sessions. Wall-clock timing on a loaded desktop drifts, so the fix is not to
%restate one against the other but to measure BOTH IN ONE INTERLEAVED ROUND,
%after which they cannot disagree.
%
%This script produces, from a single set of 300 interleaved rounds:
%   (a) every column of Table IV  -- t_it, t_fix, t_fix/t_it, M(Ts)
%   (b) every number Section III-F quotes -- dense vs factored assembly and
%       dense vs factored per-iteration product
%and it also measures the off-broadside incidence point Section III-D needs.
%
%License: GPLv2.

clear; close all; rng(2026,'twister');
addPaper2Paths();
fprintf('\n========== exp05b: Table IV and Section III-F, one round ==========\n');

p     = functionSimParams();
P_ris = functionRISPositions(p);
bin   = functionArrayResponse(P_ris, p.pBS, p.lambda);
N = p.N; K = p.K; lam = p.lambda;
gbar  = sqrt(N/K);

pU0 = zeros(3,K);
for k = 1:K
    pU0(:,k) = [p.region(1,1)+diff(p.region(1,:))*rand;
                p.region(2,1)+diff(p.region(2,:))*rand;
                p.region(3,1)+diff(p.region(3,:))*rand];
end
th0 = exp(1j*2*pi*rand(N,1));
H0  = zeros(N,K);
for k=1:K, H0(:,k) = functionUserChannel(p,P_ris,bin,pU0(:,k)); end
Cb  = functionBuildCodebook(p,P_ris,bin);

SigR  = (0.01)^2*eye(3);       % 1 cm
delta = 0.05;                  % 5 cm
vhat  = [0;1;0];               % transverse
nRep  = 300;

pf = p; pf.muTol = 0; pf.epsTol = 0; pf.Lmax = 1;
Xi0 = bsxfun(@times, conj(functionArrayResponse(P_ris,pU0,lam)), bin);
functionSteeringMatrix(H0,Cb);
pf.Imax=20;  functionBCDMM(Xi0,th0,pf);
pf.Imax=120; functionBCDMM(Xi0,th0,pf);

Fc = cell(1,K); Fp = cell(1,K); Rd = zeros(N,N,K);
for k=1:K
    Fc{k} = functionSweptMoment(P_ris,pU0(:,k),SigR,delta,vhat,lam,2,0,2);
    Fp{k} = functionRobustMomentFactor(P_ris,pU0(:,k),SigR,lam,2);
    Rd(:,:,k) = functionRobustMoment(P_ris,pU0(:,k),SigR,lam);
end

tCB=zeros(nRep,1); tDS=zeros(nRep,1); tUA=zeros(nRep,1);
tAsmD=zeros(nRep,1); tAsmP=zeros(nRep,1);
tItD=zeros(nRep,1); tItF=zeros(nRep,1); tItP=zeros(nRep,1);
tLo=zeros(nRep,1); tHi=zeros(nRep,1);
for r = 1:nRep
    tic; functionSteeringMatrix(H0,Cb); tCB(r)=toc;
    tic; Xi = bsxfun(@times,conj(functionArrayResponse(P_ris,pU0,lam)),bin); tDS(r)=toc; %#ok<NASGU>
    tic; for k=1:K, Fc{k}=functionSweptMoment(P_ris,pU0(:,k),SigR,delta,vhat,lam,2,0,2); end
    tUA(r)=toc;
    tic; for k=1:K, Rd(:,:,k)=functionRobustMoment(P_ris,pU0(:,k),SigR,lam); end
    tAsmD(r)=toc;
    tic; for k=1:K, Fp{k}=functionRobustMomentFactor(P_ris,pU0(:,k),SigR,lam,2); end
    tAsmP(r)=toc;
    tic; for k=1:K, v=Rd(:,:,k)*th0;        end, tItD(r)=toc; %#ok<NASGU>
    tic; for k=1:K, v=Fc{k}*(Fc{k}'*th0);   end, tItF(r)=toc; %#ok<NASGU>
    tic; for k=1:K, v=Fp{k}*(Fp{k}'*th0);   end, tItP(r)=toc; %#ok<NASGU>
    pf.Imax=20;  tic; functionBCDMM(Xi0,th0,pf); tLo(r)=toc;
    pf.Imax=120; tic; functionBCDMM(Xi0,th0,pf); tHi(r)=toc;
end

t_it   = median((tHi-tLo)/100);
t_set  = median(tLo) - 20*t_it;
t_itU  = t_it + median(tItF);
tfix_ds = median(tDS) + t_set;
tfix_ua = median(tDS) + median(tUA) + t_set;
phi_ds  = tfix_ds/t_it;
phi_ua  = tfix_ua/t_itU;

fprintf('\n=== TABLE IV, as measured in this round (%d interleaved rounds) ===\n', nRep);
fprintf('  Pipeline          t_it [us]   t_fix [ms]   t_fix/t_it   M at Ts=2ms\n');
fprintf('  P2-DS (nominal)  %9.2f %12.3f %12.1f %13.1f\n', ...
        t_it*1e6, tfix_ds*1e3, phi_ds, (2e-3-tfix_ds)/t_it);
fprintf('  P2-UA (moments)  %9.2f %12.3f %12.1f %13.1f\n', ...
        t_itU*1e6, tfix_ua*1e3, phi_ua, (2e-3-tfix_ua)/t_itU);
fprintf('  per-iteration penalty  %.3fx      fixed-cost penalty  %.3fx\n', ...
        t_itU/t_it, tfix_ua/tfix_ds);
fprintf('  minimum sustainable Ts: P2-DS %.3f ms   P2-UA %.3f ms\n', ...
        tfix_ds*1e3, tfix_ua*1e3);
for Ts = [1 2 3 5]*1e-3
    fprintf('    Ts = %3.0f ms -> M(DS) %8.1f   M(UA) %8.1f\n', ...
            Ts*1e3, (Ts-tfix_ds)/t_it, (Ts-tfix_ua)/t_itU);
end

Q = @(v) [median(v), prctile(v,10), prctile(v,90)]*1e3;
fprintf('\n=== SECTION III-F, same round ===\n');
q=Q(tAsmD); fprintf('  assembly  dense                    %8.3f ms  (%.3f, %.3f)\n', q);
q=Q(tUA);   fprintf('  assembly  factored (Table IV''s)    %8.3f ms  (%.3f, %.3f)\n', q);
q=Q(tAsmP); fprintf('  assembly  factored, degree 2       %8.3f ms  (%.3f, %.3f)\n', q);
q=Q(tItD);  fprintf('  product   dense                    %8.4f ms  (%.4f, %.4f)\n', q);
q=Q(tItF);  fprintf('  product   factored, compressed r=2 %8.4f ms  (%.4f, %.4f)\n', q);
q=Q(tItP);  fprintf('  product   factored, r_p = 10       %8.4f ms  (%.4f, %.4f)\n', q);
fprintf('  ---------------------------------------------------------------\n');
fprintf('  assembly  dense -> factored   %.2fx  (%.3f -> %.3f ms)\n', ...
        median(tAsmD)/median(tUA), median(tAsmD)*1e3, median(tUA)*1e3);
fprintf('  product   dense -> factored   %.1fx  (%.4f -> %.4f ms)\n', ...
        median(tItD)/median(tItF), median(tItD)*1e3, median(tItF)*1e3);
fprintf('  CONSISTENCY CHECK: t_itU - t_it = %.4f ms must equal the factored product %.4f ms\n', ...
        (t_itU-t_it)*1e3, median(tItF)*1e3);

%% off-broadside incidence at zeta = 0.08, for Section III-D
fprintf('\n=== SECTION III-D: incidence dependence at zeta = 0.08 ===\n');
zt = 0.08; nMC = 20000; rho = 25; wNF = lam*rho/p.D; sig = zt*wNF;
angs = [0 15 30 45];
lossA = zeros(numel(angs),1); ciA = zeros(numel(angs),1);
for i=1:numel(angs)
    th = angs(i)*pi/180;
    pu = rho*[cos(th); sin(th); 0];
    a0 = functionArrayResponse(P_ris,pu,lam);
    g = zeros(nMC,1);
    for m=1:nMC
        at = functionArrayResponse(P_ris, pu + sig*randn(3,1), lam);
        g(m) = abs(at'*a0)^2/N^2;
    end
    lossA(i) = 1-mean(g); ciA(i) = 1.96*std(g)/sqrt(nMC);
    fprintf('  incidence %2d deg   loss = %.4f  +- %.4f\n', angs(i), lossA(i), ciA(i));
end
cD = (p.N1+1)/(p.N1-1); pred = 2*pi^2*cD/3*zt^2;
fprintf('  closed form (27) retains its broadside value  %.4f\n', pred);
fprintf('  over-prediction: %.1f%% at broadside, %.1f%% at 45 deg\n', ...
        100*(pred/lossA(1)-1), 100*(pred/lossA(end)-1));

save('exp05b_tableIV.mat','t_it','t_itU','tfix_ds','tfix_ua','phi_ds','phi_ua', ...
     'tAsmD','tUA','tAsmP','tItD','tItF','tItP','nRep','angs','lossA','ciA','pred','p');
fprintf('\n========== done. saved exp05b_tableIV.mat ==========\n');
