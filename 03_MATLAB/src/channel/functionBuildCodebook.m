function [C, Qcb] = functionBuildCodebook(p, P_ris, bin)
%Build the region-gridded near-field codebook C_NF of spherical-wavefront
%steering (focusing) vectors, as in Sec. III-B of the paper.
%
%Each codeword focuses the RIS at a grid point q:
%   c_s = conj(a(q_s)) .* b_in      (unit-modulus per element)
%so that |h_k^H c_s| is maximised when q_s coincides with user k.
%
%INPUT:
% p      = parameter struct
% P_ris  = 3 x N RIS positions
% bin    = N x 1 incident steering vector b_in = a(p_BS)
%
%OUTPUT:
% C      = N x S codebook matrix (S = cbNx*cbNy*cbNz)
% Qcb    = 3 x S focal-point coordinates [m]
%
%This is version 1.0 (Last edited: 2026-06-02)
%License: GPLv2.

xg = linspace(p.region(1,1), p.region(1,2), p.cbNx);
yg = linspace(p.region(2,1), p.region(2,2), p.cbNy);
zg = linspace(p.region(3,1), p.region(3,2), p.cbNz);
[X,Y,Z] = ndgrid(xg, yg, zg);
Qcb = [X(:).'; Y(:).'; Z(:).'];            %3 x S

A = functionArrayResponse(P_ris, Qcb, p.lambda);   %N x S, a(q_s)
C = conj(A) .* bin;                                %focus + incident
end
