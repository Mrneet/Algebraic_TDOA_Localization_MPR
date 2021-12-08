function [ bias ] = Bias_GTRS_MPR( senPos, srcLoc, Q )
%function [ bias ] = Bias_GTRS_MPR( senPos, srcLoc, Q )
% Evaluating the theoretical bias of the source location estimate in MPR 
% for the GTRS method.  It can work with 2-D or 3-D scanario.
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
g0 = 1/r(1);
R = diag(rd);

B1 = -diag(r(2:end)/r(1));
W = eye(M-1)/(B1*Q*B1);

q = diag(Q);
Go = [senPos(:,2:end)', 0.5*(rd.^2 - sum(senPos(:,2:end)'.^2, 2))];
phi = [srcLoc*g0;g0];
S = diag([ones(N,1);0]);
iP1o = eye(N+1)/(Go'*W*Go);
id1o = 1/(phi'*S*iP1o*S*phi);
H3 = S*(phi*phi')*S;

P6 = iP1o*H3*iP1o;
H4 = iP1o - id1o*P6;
bias1 = -g0/2*H4*Go'*W*q + H4*[zeros(N,1);trace(R*W*B1*Q)];
P7 = W*Go*iP1o*Go'*W;
H1 = iP1o*Go'*W;
bias2 = -H4*[zeros(N,1);trace(R*P7*B1*Q)] - H4*Go'*W*R*Q*B1*H1(end,:)';
P8 = W*Go*P6*Go'*W;
H5 = P6*Go'*W;
bias3 = id1o*iP1o*[zeros(N,1);trace(R*P8*B1*Q)] + id1o*iP1o*Go'*W*R*Q*B1*H5(end,:)';
p1 = iP1o*S*phi;
bias4 = -2*id1o^2*p1(end)*P6*Go'*W*B1*Q*R*W*Go*p1;
Bias_p = bias1 + bias2 + bias3 + bias4;

if N == 2 % 2-D
    doa = atan2(srcLoc(2),srcLoc(1));
    D1 = [-sin(doa),cos(doa),0;
          0,        0,       1];
elseif N == 3 % 3-D
    theta = atan2(srcLoc(2), srcLoc(1));
    phi = atan2(srcLoc(3), norm(srcLoc(1:2),'fro'));
    D1 = [-sin(theta)/cos(phi), cos(theta)/cos(phi),    0,          0;
          -cos(theta)*sin(phi), -sin(theta)*sin(phi),   cos(phi),   0;
          0,                    0,                      0,          1];
else
    error('Please check your input format of sensor positions');
end

bias = D1*Bias_p;