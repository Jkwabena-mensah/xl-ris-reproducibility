function h = functionUserChannel(p, P_ris, bin, pk)
%Effective near-field reflecting channel h_k from the BS to a user via the
%XL-RIS, now with a LoS + NLOS (near-field Rician) reflected steering vector
%[eq.(1)-(2) of the paper]:
%
%   b_re,k = sqrt(kappa/(kappa+1)) * a(p_k)
%          + sqrt(1/(kappa+1)) * (1/sqrt(L)) * sum_l beta_{k,l} a(q_{k,l})
%   h_k    = alpha * ( conj(b_re,k) .* b_in ),   alpha = 1/N
%
%where a(.) is the spherical-wavefront response, q_{k,l} are L near-field
%scatterers around the user, and beta_{k,l} ~ CN(0,1).  As kappa -> inf this
%reduces to the pure-LoS spherical-wavefront model.
%
%The direct BS->user link is treated separately (functionDirectChannel.m);
%for the blocked-obstacle case (p.rhoDirect = 0) it is zero and does not
%enter the rate, so this function alone suffices for the simulations.
%
%INPUT:
% p      = parameter struct (uses lambda, N, alpha, kappa, Lscat, scatSpread)
% P_ris  = 3 x N RIS positions
% bin    = N x 1 incident steering vector
% pk     = 3 x 1 user position [m]
%
%OUTPUT:
% h      = N x 1 effective reflected channel vector
%
%This is version 2.0 (Last edited: 2026-06-02)  -- adds NLOS multipath
%License: GPLv2.

aLoS = functionArrayResponse(P_ris, pk, p.lambda);     %N x 1 LoS response

L  = p.Lscat;
bN = zeros(p.N,1);
for l = 1:L
    %near-field scatterer in the vicinity of the user (smaller z-spread)
    q    = pk + p.scatSpread*[randn; randn; 0.3*randn];
    beta = (randn + 1i*randn)/sqrt(2);                 %CN(0,1) gain
    bN   = bN + beta*functionArrayResponse(P_ris, q, p.lambda);
end

bre = sqrt(p.kappa/(p.kappa+1))*aLoS ...
    + sqrt(1/(p.kappa+1))*(1/sqrt(L))*bN;              %Rician steering

h = p.alpha * (conj(bre) .* bin);
end
