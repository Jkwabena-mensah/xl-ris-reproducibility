function A = functionArrayResponse(P_ris, Q, lambda)
%Near-field (spherical-wavefront) array-response vectors from the RIS
%elements to one or more spatial points.
%
%The n-th entry for point q is  exp(-j*2*pi*||p_n - q|| / lambda).
%
%INPUT:
% P_ris  = 3 x N RIS element positions [m]
% Q      = 3 x S matrix of spatial points (users or focal points) [m]
% lambda = wavelength [m]
%
%OUTPUT:
% A      = N x S complex matrix; column s is the response to point Q(:,s)
%
%This is version 1.0 (Last edited: 2026-06-02)
%License: GPLv2.

N = size(P_ris,2);
S = size(Q,2);
A = zeros(N,S);
for s = 1:S
    d = sqrt(sum((P_ris - Q(:,s)).^2, 1));   %1 x N distances
    A(:,s) = exp(-1j*2*pi*d(:)/lambda);
end
end
