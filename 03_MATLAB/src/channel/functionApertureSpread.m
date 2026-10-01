function [C, ubar, U] = functionApertureSpread(P_ris, pUser)
%Compute the aperture spread matrix C of the near-field array response at a
%given user location [eq. (11) of the article].
%
%The unit vector from RIS element n toward the user is
%   u_n = (p - p_n^RIS) / || p - p_n^RIS ||,
%and with ubar the aperture average of u_n and utilde_n = u_n - ubar,
%   C = (1/N) * sum_n utilde_n * utilde_n'.
%
%C is the covariance of the element-to-user directions across the aperture. It
%governs the beam-focusing loss caused by a position error (Theorem 1), whereas
%C + ubar*ubar' governs the change in the desired response and hence the work
%demanded of the optimiser (Proposition 2). The difference between the two is
%exactly the rank-one term ubar*ubar', which is why a radial position error
%costs effort but almost no focusing gain.
%
%In the far field all u_n coincide, so C -> 0 and a position error becomes a
%phase common to every element, costing no array gain at all.
%
%INPUT:
% P_ris = 3 x N matrix of RIS element positions [m]
% pUser = 3 x 1 user position [m]
%
%OUTPUT:
% C     = 3 x 3 aperture spread matrix (dimensionless, symmetric PSD)
% ubar  = 3 x 1 mean element-to-user direction (unit-norm to leading order)
% U     = N x 3 matrix whose n-th row is u_n'
%
%This is version 1.0 (Last edited: 2026-09-07)
%License: GPLv2.

N = size(P_ris,2);

V = pUser - P_ris;                              %3 x N, element -> user vectors
U = (V ./ vecnorm(V,2,1)).';                    %N x 3, unit directions (rows)

ubar  = mean(U,1).';                            %3 x 1 mean direction
Utild = U - ubar.';                             %N x 3 centred directions

C = (Utild.' * Utild) / N;                      %3 x 3
C = (C + C.')/2;                                %Symmetrise against round-off
end
