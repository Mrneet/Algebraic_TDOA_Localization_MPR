function [ CRLB ] = CRLB_MPR( senPos, srcLoc, Q )
%function [ CRLB ] = CRLB_MPR( senPos, srcLoc, Q )
% Evaluating the CRLB for the source location in MPR.  It can work with 
% 2-D or 3-D scanario.
%
% Input:
%   senPos:     NxM, positions of reciving sensors, each column is a sensor 
%               position and the first column is the reference sensor
%               position;
%   srcLoc:     Mx1, source location 
%   Q:          (M-1)x(M-1), the covariance matrix of range differences (TDOAs).
% Output:
%   CRLB:       NxN, CRLB for source location in MPR
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
u0 = srcLoc/r(1);

p = repmat(u0,1,M-1) - senPos(:,2:end)/r(1);
P = p./repmat(sqrt(sum(p.^2,1)),size(p,1),1);

if N == 2
    theta = atan2(srcLoc(2),srcLoc(1));
    MM = r(1)* ...
        [-sin(theta),  -cos(theta)*r(1);
         cos(theta),   -sin(theta)*r(1)];
elseif N == 3
    theta = atan2(srcLoc(2), srcLoc(1));
    phi = atan2(srcLoc(3), norm(srcLoc(1:2),'fro'));
    MM = r(1)* ...
        [-sin(theta)*cos(phi),  -cos(theta)*sin(phi),   -cos(theta)*cos(phi)*r(1);
         cos(theta)*cos(phi),   -sin(theta)*sin(phi),   -sin(theta)*cos(phi)*r(1);
         0,                     cos(phi),               -sin(phi)*r(1)];
else
    error('Please check your input format of sensor positions');
end

L = [zeros(M-1,N-1),ones(M-1,1)*r(1)^2];
U = P'*MM + L;
CRLB = inv(U'/Q*U);

