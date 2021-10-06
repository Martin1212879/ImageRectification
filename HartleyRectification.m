function [H1,H2] = HartleyRectification(F,mu,nu,He2)
%% Based on Theory and Practice of Projective Rectification-RICHARD I. HARTLEY
%
% HartleyRectfication return H1 & H2 homography matrix 
%
% Syntaxis: function [H1,H2] = HartleyRectification(F,mu,nu,He2)
%
% Input:
%       F - Fundamentala Matrix (3x3)
%       mu,nu - matching points image 1 & 2 respectively (2xn)
%       He2 - Epipolar point image 2 (3x1)
% 
% Output:
%       H1 & H2 homography matrix image 1 & 2 respectively (3x3)

%% Rotation
He2 = He2/He2(end);
c = He2(1)/sqrt(He2(1)^2 + He2(2)^2);
s = He2(2)/sqrt(He2(1)^2 + He2(2)^2);
R = [c   s  0;
     -s  c  0;
     0   0  1];
p = R*He2;
%% Epipole to infinity  
G = [1          0   0;
     0          1   0;
     -1/p(1)    0   1]; 
H2 = G*R;  
%% Matching transformation 
v = [1;1;1];
M = mcross(He2)*F + He2*v';
H0 = H2*M;
%% affine transformation 
uig = H0*hom(mu);
uigp = H2*hom(nu);
 
% Linear least-squares 
AA = uig';
YY = uigp(1,:)';
theta = inv(AA'*AA)*AA'*YY;

a1 = theta(1);
a2 = theta(2);
a3 = theta(3);

A = [a1 a2 a3;
      0  1  0;
      0  0  1];
H1 = A*H0;
