function P = functionRISPositions(p)
%Compute the 3-D positions of all XL-RIS elements.
%The RIS is a uniform planar array (UPA) on the yOz-plane, centred at the
%origin, with inter-element spacing p.dElem.
%
%INPUT:
% p = parameter struct from functionSimParams()
%
%OUTPUT:
% P = 3 x N matrix; column n is the [x;y;z] position of element n [m]
%
%This is version 1.0 (Last edited: 2026-06-02)
%License: GPLv2.

yIdx = (0:p.N2-1) - (p.N2-1)/2;     %Centred column indices
zIdx = (0:p.N1-1) - (p.N1-1)/2;     %Centred row indices

[Y, Z] = meshgrid(yIdx*p.dElem, zIdx*p.dElem);   %N1 x N2 grids
P = [zeros(1, p.N); Y(:).'; Z(:).'];             %3 x N (x = 0 plane)
end
