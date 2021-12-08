function [ bias ] = Bias_SUM_MPR( senPos, srcLoc, Q )
%function [ bias ] = Bias_SUM_MPR( senPos, srcLoc, Q )
% Evaluating the theoretical bias of the source location estimate in MPR 
% for the two-stage WLS method.  It can work with 2-D or 3-D scanario.
%
% Input:
%   senPos:     NxM, positions of reciving sensors, each column is a sensor 
%               position and the first column is the reference sensor
%               position;
%   srcLoc:     Mx1, source location 
%   Q:          (M-1)x(M-1), the covariance matrix of range differences (TDOAs).
% Output:
%   bias:       Nx1, bias in [theta; phi(absent in 2-D); g]
%
% Reference: Y. Sun, K. C. Ho, and Q. Wan, "Solution and analysis of TDOA 
%  localization of a near or distant source in closed-form," IEEE Trans. 
%  Signal Process., vol. 67, no. 2, pp. 320-335, Jan. 2019.
%
% Yimao Sun, K. C. Ho    02-28-2019
%
%       Copyright (C) 2019
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

[N,M] = size(senPos);
r = sqrt(sum((senPos-repmat(srcLoc,1,M)).^2,1))';
rd = r(2:end)-r(1);
R = diag(rd);
u0 = srcLoc/r(1);

B1o = -diag(r(2:end)/r(1));
W1 = inv(B1o*Q*B1o);

G1o = [senPos(:,2:end)', 0.5*(rd.^2 - sum(senPos(:,2:end)'.^2, 2))];
P1o = G1o'*W1*G1o;
C1 = inv(P1o);    % convariance of first stage
H1 = C1*G1o'*W1;

E10 = -0.5/r(1)*H1*diag(Q);
E11 = C1*[zeros(N,1);trace(R*W1*B1o*Q)];
E12 = -H1*R*Q*B1o*H1(N+1,:)';
E13 = -C1*[zeros(N,1);trace(R*W1*G1o*H1*B1o*Q)];
bias1 = E10 + E11 + E12 + E13;

B2o = diag([2*u0;1]);
iB2o = inv(B2o);
W2o = iB2o*(G1o'*W1*G1o)*iB2o;
G2 = [eye(N-1),zeros(N-1,1);-ones(1,N-1),0;zeros(1,N-1),1];
P2o = G2'*W2o*G2;
H2 = (G2'*W2o*G2)\G2'*W2o;

T = diag([ones(1,N),0]);
E20 = H2*T*diag(C1);
E21 = H2*B2o*bias1;
P3 = eye(N+1) - G2*H2;
P4 = iB2o*P3*B2o*H1*B1o;
alp = iB2o*G1o'*W1*R*Q*P4(N+1,:)' + iB2o*[zeros(N,1);trace(R*W1*G1o*P4*Q)];
P5 = P3*B2o*C1;
P5W = W2o*P5;
bet = -2*iB2o*T*diag(P5W);
gam = -2*W2o*iB2o*T*diag(P5);
E22 = P2o\G2'*(alp + bet + gam);
bias2 = E20 + E21 + E22;

if N == 2 % 2-D
    doa = atan2(srcLoc(2),srcLoc(1));
    D2 = [-sin(2*doa),0;0,1];
elseif N == 3 % 3-D
    theta = atan2(srcLoc(2), srcLoc(1));
    phi = atan2(srcLoc(3), norm(srcLoc(1:2),'fro'));
    D2 = [-2*sin(theta)*cos(theta)*cos(phi)^2,    -2*cos(theta)^2*sin(phi)*cos(phi),  0;
          2*sin(theta)*cos(theta)*cos(phi)^2,     -2*sin(theta)^2*sin(phi)*cos(phi),	0;
          0,                                      0,                                  1];
else
    error('Please check your input format of sensor positions');
end

bias = D2\bias2;
