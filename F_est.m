
function [F,He1,He2] = F_est(mu,nu)
% F_est return 3x3 Fundamental matrix 
% 
% Syntaxys: [F,He1,He2] = F_est(mu,nu)
% 
% Input:
%        mu,nu are matching points image 1 & 2 respectively (2xn)
%
% Output: 
%        F is a fundamental matrix (3x3)
%        He1, He2 are epipoles in homogeneous coordinates

%% Size of points
N = size(mu,2);
%% Matching points  (nx1)
x1 = mu(1, :)';
y1 = mu(2, :)';
x2 = nu(1, :)';
y2 = nu(2, :)';
%% Regresion Matrix (nx8)
A = [x2 .* x1, x2 .* y1, x2, y2 .* x1, y2 .* y1, y2, x1, y1 ones(N,1)];
%% SVD
[~, ~, V] = svd(A);
fMatrix = [V(1, 9), V(2, 9), V(3, 9); V(4, 9), V(5, 9), V(6, 9); V(7, 9), V(8, 9), V(9, 9)];
%% Rank 2
[U, S, V] = svd(fMatrix);
S(end)=0;
F = U*S*V';
%% Epipoles
[He1,He2] = epipole_est(F);
