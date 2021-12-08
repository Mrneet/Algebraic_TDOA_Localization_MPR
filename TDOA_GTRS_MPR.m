function [ varargout ] = TDOA_GTRS_MPR( senPos, rd, Q )
%function [ varargout ] = TDOA_GTRS_MPR( senPos, rd, Q )
% Estimating the source location in MPR (angle, inverse-range) by the 
% closed-form solution using GTRS.  It can work with 2-D or 3-D 
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

[N,M]=size(senPos);     % localization dimension, number of sensors

h1 = -rd;
l = sqrt(sum((senPos(:,2:end) - repmat(senPos(:,1),1,M-1)).^2, 1))';
s_bar = senPos(:,2:end)./repmat(l',N,1);

% first stage
G1 = [senPos(:,2:end)', 0.5*(rd.^2 - sum(senPos(:,2:end)'.^2, 2))];
B1 = eye(M-1);
W = eye(M-1)/(B1*Q*B1);

for iter = 1:2,     % repeat once after weighting matrix update
    A = G1(:,1:N);
    a = G1(:,N+1);
    O = W-W*a/(a'*W*a)*a'*W;
    AOA = A'*O*A;
    S = eye(N);
    
    [U,D,~] = svd(AOA);
    r = diag(D);
    k = U'*A'*O*h1;
    
    if N == 2
        x(1) = 1;
        x(2) = 2*r(1)+2*r(2);
        x(3) = r(1)^2+2*r(1)*r(2)+r(2)^2-k(1)^2-k(2)^2;
        x(4) = 2*r(1)^2*r(2)+2*r(1)*r(2)^2-2*k(1)^2*r(2)-2*k(2)^2*r(1);
        x(5) = r(1)^2*r(2)^2-k(1)^2*r(2)^2-k(2)^2*r(1)^2;
    elseif N == 3
        x(1) = 1;
        x(2) = 2*r(3)+2*r(2)+2*r(1);
        x(3) = r(3)^2+4*r(2)*r(3)+4*r(1)*r(3)+r(2)^2+4*r(1)*r(2)+r(1)^2-k(3)^2-k(2)^2-k(1)^2;
        x(4) = 2*r(2)*r(3)^2+2*r(1)*r(3)^2+2*r(2)^2*r(3)+8*r(1)*r(2)*r(3)+2*r(1)^2*r(3)-2*k(2)^2*r(3)-2*k(1)^2*r(3)+2*r(1)*r(2)^2+2*r(1)^2*r(2)-2*k(3)^2*r(2)-2*k(1)^2*r(2)-2*k(3)^2*r(1)-2*k(2)^2*r(1);
        x(5) = r(2)^2*r(3)^2+4*r(1)*r(2)*r(3)^2+r(1)^2*r(3)^2-k(2)^2*r(3)^2-k(1)^2*r(3)^2+4*r(1)*r(2)^2*r(3)+4*r(1)^2*r(2)*r(3)-4*k(1)^2*r(2)*r(3)-4*k(2)^2*r(1)*r(3)+r(1)^2*r(2)^2-k(3)^2*r(2)^2-k(1)^2*r(2)^2-4*k(3)^2*r(1)*r(2)-k(3)^2*r(1)^2-k(2)^2*r(1)^2;
        x(6) = 2*r(1)*r(2)^2*r(3)^2+2*r(1)^2*r(2)*r(3)^2-2*k(1)^2*r(2)*r(3)^2-2*k(2)^2*r(1)*r(3)^2+2*r(1)^2*r(2)^2*r(3)-2*k(1)^2*r(2)^2*r(3)-2*k(2)^2*r(1)^2*r(3)-2*k(3)^2*r(1)*r(2)^2-2*k(3)^2*r(1)^2*r(2);
        x(7) = r(1)^2*r(2)^2*r(3)^2-k(1)^2*r(2)^2*r(3)^2-k(2)^2*r(1)^2*r(3)^2-k(3)^2*r(1)^2*r(2)^2;
    end
    root = roots(x);
    
    % delete complex roots
    reRoot = root(imag(root)==0);
    L = length(reRoot);
    if L == 0 % to insure that Y is not empty
        [~,I] = min(imag(root));
        reRoot = real(root(I));
        L = 1;
    end
    
    Y = zeros(N+1,L);J = zeros(1,L);
    for i = 1:L
        Y(1:N,i) = (AOA+reRoot(i)*S)\A'*O*h1;
        Y(N+1,i) = (a'*W*a)\a'*W*(h1-A*Y(1:N,i));
        J(i) = (h1-G1*Y(:,i))'*W*(h1-G1*Y(:,i));
    end
    [~,ind] = min(J);
    Phi0 = Y(:,ind);
    Phi = sign(real(Phi0)).*abs(Phi0);   % to keep it real
    
    b = 1 + l.^2*Phi(end)^2 - 2*Phi(end)*l.*(s_bar'*Phi(1:N));
    B1 = -diag(sqrt(b));
    W = eye(M-1)/(B1*Q*B1);
end

if N == 2
    theta = atan2(Phi(2),Phi(1));
    g = Phi(3);
    pos = [cos(theta);sin(theta)]/g;
    varargout{1} = theta;
    varargout{2} = g;
    varargout{3} = pos;
elseif N == 3
    theta = atan2(Phi(2),Phi(1));
    phi = atan2(Phi(3),sqrt(Phi(1)^2+Phi(2)^2));
    g = Phi(4);
    pos = [cos(theta)*cos(phi);sin(theta)*cos(phi);sin(phi)]/g;
    varargout{1} = theta;
    varargout{2} = phi;
    varargout{3} = g;
    varargout{4} = pos;
else
    error('Please check your input format of sensor positions');
end


