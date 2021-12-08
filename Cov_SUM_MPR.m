function [ cov2WLS ] = Cov_SUM_MPR( senPos, srcLoc, Q )
%function [ cov2WLS ] = Cov_SUM_MPR( senPos, srcLoc, Q )
% Evaluating the theoretical covariance matrix of the source location estimate
% in MPR for the two-stage WLS method.  It can work with 2-D or 3-D scanario.
%
% Input:
%   senPos:     NxM, positions of reciving sensors, each column is a sensor 
%               position and the first column is the reference sensor
%               position;
%   srcLoc:     Mx1, source location 
%   Q:          (M-1)x(M-1), the covariance matrix of range differences (TDOAs).
% Output:
%   cov2WLS:    NxN, covariance of source location estimate in MPR
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

l = sqrt(sum((senPos(:,2:end) - repmat(senPos(:,1),1,M-1)).^2, 1))';
s_bar = senPos(:,2:end)./repmat(l',N,1);
u0 = srcLoc/r(1);
b = 1 + l.^2/r(1)^2 - 2/r(1)*l.*(s_bar'*u0);
B1 = -diag(sqrt(b));

G1 = [senPos(:,2:end)', 0.5*(rd.^2 - sum(senPos(:,2:end)'.^2, 2))];

C1 = inv(G1'/(B1*Q*B1)*G1);    % convariance of first stage

B2 = diag([2*u0;1]);
G2 = [eye(N-1),zeros(N-1,1);-ones(1,N-1),0;zeros(1,N-1),1];
C2 = inv(G2'/(B2*C1*B2)*G2);

if N == 2 % 2-D
    doa = atan2(srcLoc(2),srcLoc(1));
    D2 = [-sin(2*doa),0;0,1];
elseif N == 3 % 3-D
    theta = atan2(srcLoc(2), srcLoc(1));
    phi = atan2(srcLoc(3), norm(srcLoc(1:2),'fro'));
    D2 = [-2*sin(theta)*cos(theta)*cos(phi)^2,    -2*cos(theta)^2*sin(phi)*cos(phi),  0;
          2*sin(theta)*cos(theta)*cos(phi)^2,     -2*sin(theta)^2*sin(phi)*cos(phi),  0;
          0,                                      0,                                  1];
else
    error('Please check your input format of sensor positions');
end

cov2WLS = D2\C2/D2;