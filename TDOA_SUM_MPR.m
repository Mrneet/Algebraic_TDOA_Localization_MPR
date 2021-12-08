function [ varargout ] = TDOA_SUM_MPR( senPos, rd, Q )
%function [ varargout ] = TDOA_SUM_MPR( senPos, rd, Q )
% Estimating the source location in MPR (angle, inverse-range) by the 
% closed-form solution using two-stage WLS.  It can work with 2-D or 3-D 
% scanario.
%
% Input:
%   senPos:     NxM, positions of reciving sensors, each column is a sensor 
%               position and the first column is the reference sensor
%               position;
%   rd:         (M-1)x1, range difference (TDOA) measurement vector;
%   Q:          (M-1)x(M-1), covariance matrix of range differences (TDOAs).
% Output:
%   varargout: including
%       theta:	azimuth estimate
%       phi:	elevation estimate (absent for 2-D)
%       g:      g estimate
%       pos:	Nx1, source position estimate
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

FirstStageOnly=0;       % set to one if we do not apply second stage
[N,M]=size(senPos);     % localization dimension, number of sensors

h1 = -rd;
l = sqrt(sum((senPos(:,2:end) - repmat(senPos(:,1),1,M-1)).^2, 1))';
s_bar = senPos(:,2:end)./repmat(l',N,1);

% -- first stage --
G1 = [senPos(:,2:end)', 0.5*(rd.^2 - sum(senPos(:,2:end)'.^2, 2))];
Phi1 = (G1'/Q*G1)\(G1')/Q*h1;   % initial estimation

b = 1 + l.^2*Phi1(end)^2 - 2*Phi1(end)*l.*(s_bar'*Phi1(1:N));
B1 = -diag(sqrt(b));
W1_inv = B1*Q*B1;
Phi20 = ((G1'/W1_inv*G1)\G1'/W1_inv*h1); % solution of 1st stage
Phi2 = sign(real(Phi20)).*abs(Phi20);    % to keep it real

% results from first stage
Psi=Phi2;

if (~FirstStageOnly)
    % -- second stage --
    B2 = diag([2*Phi2(1:N);1]);
    h2 = [Phi2(1:N-1).^2;Phi2(N)^2-1;Phi2(end)];
    G2 = [eye(N-1),zeros(N-1,1);-ones(1,N-1),0;zeros(1,N-1),1];
    W2 = B2\(G1'/W1_inv*G1)/B2;
    Psi0 = ((G2'*W2*G2)\(G2')*W2*h2);   % solution of 2nd stage
    Psi = sign(real(Psi0)).*abs(Psi0);  % to keep it real
end;

if N == 2
    theta = atan2( abs(sqrt(1-Psi(1)))*sign(Phi2(2)), abs(sqrt(Psi(1)))*sign(Phi2(1)) );
    g = Psi(2);
    pos = [cos(theta);sin(theta)]/g;
    varargout{1} = theta;
    varargout{2} = g;
    varargout{3} = pos;
elseif N == 3
    theta = atan2( abs(sqrt(Psi(2)))*sign(Phi2(2)), abs(sqrt(Psi(1)))*sign(Phi2(1)) );
    phi = acos(sqrt(sum(Psi(1:2)))*sign(Phi2(3)));
    g = Psi(end);
    pos = [cos(theta)*cos(phi);sin(theta)*cos(phi);sin(phi)]/g;
    
    if (FirstStageOnly)
        theta = atan2(Phi2(2),Phi2(1));   % result of 1st stage
        phi = atan2(Phi2(3),sqrt(Phi2(1)^2+Phi2(2)^2));
        g = Phi2(4);
        pos = [cos(theta(1))*cos(phi(1));sin(theta(1))*cos(phi(1));sin(phi(1))]/g(1);
    end; 
    varargout{1} = theta;
    varargout{2} = phi;
    varargout{3} = g;
    varargout{4} = pos;
else
    error('Please check your input format of sensor positions');
end

