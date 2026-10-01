function rho = functionFocusingGain(P_ris, pDesign, e, aDesign, lambda)
%Normalised beam-focusing gain obtained when the RIS is configured to focus at
%pDesign while the user actually occupies pDesign + e [eq. (20) of the article]:
%
%   rho = | a(pDesign + e)^H * a(pDesign) |^2 / N^2,   rho in [0,1].
%
%The array responses are evaluated EXACTLY from the spherical-wave model, with
%no expansion, so this function provides the ground truth against which the
%second-order predictions of Theorem 1 and Corollary 1 are tested.
%
%INPUT:
% P_ris   = 3 x N RIS element positions [m]
% pDesign = 3 x 1 position the surface is configured for [m]
% e       = 3 x 1 position error, so the true user position is pDesign + e [m]
% aDesign = N x 1 precomputed response a(pDesign); pass [] to compute it here
% lambda  = wavelength [m]
%
%OUTPUT:
% rho     = scalar normalised focusing gain in [0,1]
%
%This is version 1.0 (Last edited: 2026-09-07)
%License: GPLv2.

N = size(P_ris,2);

if isempty(aDesign)
    aDesign = functionArrayResponse(P_ris, pDesign, lambda);
end

aTrue = functionArrayResponse(P_ris, pDesign + e, lambda);   %N x 1

rho = abs(aTrue' * aDesign)^2 / N^2;
end
