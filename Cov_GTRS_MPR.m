function [ covGTRS ] = Cov_GTRS_MPR( senPos, srcLoc, Q )
%function [ covGTRS ] = Cov_GTRS_MPR( senPos, srcLoc, Q )
% Evaluating the theoretical covariance matrix of the source location 
% estimate in MPR for the GTRS method.  It can work with 2-D or 3-D scanario.
%
% Input:
%   senPos:     NxM, positions of reciving sensors, each column is a sensor 
%               position and the first column is the reference sensor
%               position;
%   srcLoc:     Mx1, source location 
%   Q:          (M-1)x(M-1), the covariance matrix of range differences (TDOAs).
% Output:
%   covGTRS:    NxN, covariance of source location estimate in MPR
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
W = eye(M-1)/(B1*Q*B1);

G1 = [senPos(:,2:end)', 0.5*(rd.^2 - sum(senPos(:,2:end)'.^2, 2))];

phi = [u0;1/r(1)];
C = inv(G1'*W*G1);
S = diag([ones(N,1);0]);
Cov2 = C*(eye(N+1)-(S*phi*phi'*S*C)/(phi'*S*C*S*phi));

if N == 2 % 2-D
    doa = atan2(srcLoc(2),srcLoc(1));
    D1 = [-sin(doa),cos(doa),0;
          0,0,1];
elseif N == 3 % 3-D
    theta = atan2(srcLoc(2), srcLoc(1));
    phi = atan2(srcLoc(3), norm(srcLoc(1:2),'fro'));
    D1 = [-sin(theta)/cos(phi), cos(theta)/cos(phi),    0,          0;
          -cos(theta)*sin(phi), -sin(theta)*sin(phi),   cos(phi),   0;
          0,                    0,                      0,          1];
else
    error('Please check your input format of sensor positions');
end

covGTRS = D1*Cov2*D1';