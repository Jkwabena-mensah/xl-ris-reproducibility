function Xi = functionSteeringMatrix(H, C)
%Select, for each user, the best near-field codeword (codebook search),
%and assemble the steering matrix Xi = [xi_1,...,xi_K]  (Sec. III-B).
%
%   xi_k = argmax_{c in C} |h_k^H c|^2
%
%INPUT:
% H = N x K matrix of user channels (column k is h_k)
% C = N x S codebook matrix
%
%OUTPUT:
% Xi = N x K steering matrix
%
%This is version 1.0 (Last edited: 2026-06-02)
%License: GPLv2.

K = size(H,2);
N = size(H,1);
Xi = zeros(N,K);
metric = abs(H' * C).^2;          %K x S correlation magnitudes
[~, idx] = max(metric, [], 2);    %best codeword index per user
for k = 1:K
    Xi(:,k) = C(:, idx(k));
end
end
