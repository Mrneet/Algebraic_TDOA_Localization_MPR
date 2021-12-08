% This program generates the localization result of a source in MPR with
% TDOA measurements by the algebraic closed-form solutions (2WLS and GTRS)
%
% The figures generated from this code correspond to
%   Figs. 2-3 (Accuracy_vs_Range=1), or
%   Figs. 4-7 (Accuracy_vs_Range=0)
% in the reference paper.
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

clear all; 
%close all; clc;

warning off;
clor = [0, 114, 189;217, 83, 25;237, 177, 32;126, 47, 142;119, 172, 48;77, 190, 238;162, 20, 47]/256;
rng('default');

% *****************************************************************
Accuracy_vs_Range=1;      % 1 for Figs. 2-3, 0 for Figs. 4-7
Case3D=1;                 % 1 if 3D, 0 if 2D
% *****************************************************************


% ----- simulation setting -----
% sensor position
senPos = [
    0         0         0
    0.4664   -0.8712    0.6294
    1.2402   -0.0798    0.8116
    0.4437    0.6323   -0.7460
    0.5798   -0.2362    0.8268
    1.0502   -0.5172    0.2647
    -0.8156    1.0448   -0.8049]' * 1e3;

if (~Case3D)
    senPos=senPos(1:2,:);
end;

% target direction
theta = 22.13*pi/180;
phi = 14.41*pi/180;

[N,M] = size(senPos);

if (N==2) phi=0; end;

if (Accuracy_vs_Range)
    sigmaSquareDB = 0;  % 10log(m^2)
    range = [5,10,(20:20:460)]*1e3;  % m
    NumEnsembles = 2000;
else
    sigmaSquareDB = -40:10:50;  % 10log(m^2)
    range = 15*1e3;         % m
    NumEnsembles = 1000;
end;

K =  length(sigmaSquareDB);     % number of noise levels
R = length(range);              % number of ranges

if (N==3)
    srcLoc=13e3*[cos(theta)*cos(phi); sin(theta)*cos(phi); sin(phi)]';
    figure(100); plot3(senPos(1,:),senPos(2,:),senPos(3,:),'o','markersize',8,'linewidth',1);
    hold on;xlabel('x (km)','fontsize',12);ylabel('y (km)','fontsize',12);;zlabel('z (km)','fontsize',12);
    plot3(srcLoc(1), srcLoc(2), srcLoc(3), 'x','markersize',8,'linewidth',1);grid on;hold off;
    title('Localization Geometry, Sensors(circle), Source(cross)');
else
    srcLoc=13e3*[cos(theta); sin(theta)]';
    figure(100); plot(senPos(1,:),senPos(2,:),'o','markersize',8,'linewidth',1);
    hold on;xlabel('x (km)','fontsize',12);ylabel('y (km)','fontsize',12);
    plot(srcLoc(1), srcLoc(2), 'x','markersize',8,'linewidth',1);grid on;hold off;
    title('Localization Geometry, Sensors(circle), Source(cross)');
end;

% ----- Monte-Carlo Simulation -----
[eTh1,eTh2,eTh3,eTh4,eTh5,eTh6, ePh1,ePh2,ePh3,ePh4,ePh5,ePh6, eg1,eg2,eg3,eg4,eg5,eg6, ...
    ep1,ep2,ep3,ep4,ep5,ep6, er1,er2,er3,er4,er5,er6, uTh1,uTh2,uTh3,uTh4,uTh5,uTh6, ...
    uPh1,uPh2,uPh3,uPh4,uPh5,uPh6, ug1,ug2,ug3,ug4,ug5,ug6] = deal(zeros(K,NumEnsembles));

aveNse = 0;
for l=1:NumEnsembles
    aveNse = aveNse + randn(M,1);
end
aveNse = aveNse/NumEnsembles/sqrt(2);
PP = aveNse(2:end) - aveNse(1);

disp('Simulation is running ...');

for ir = 1:R,   % loop through ranges
    disp(['Range: ',num2str(range(ir)/1000),'km, ',num2str(ir),'/',num2str(R),' ...']);
    % -- Generate Data --
    % source location
    if (N==3)
        srcLoc = range(ir) * [cos(theta)*cos(phi); sin(theta)*cos(phi); sin(phi)];
    else
        srcLoc = range(ir) * [cos(theta); sin(theta)];
    end;
    % true range
    r = sqrt(sum((repmat(srcLoc,1,M)-senPos).^2,1))';
    % true TDOAs
    rd = r(2:end) - r(1);
    
    for k = 1:K,    % loop through noise powers
        disp(['Noise power (10log(\sigma^2): ',num2str(sigmaSquareDB(k)),', ',num2str(k),'/',num2str(K),' ...']);
        Q = 10^(sigmaSquareDB(k)/10) * (ones(M-1, M-1)+eye(M-1))/2;
        
        % Calculate CRLB
        CRB = TDOALocCRLB_MPR( senPos, srcLoc, Q );
        if (N==3)
            CRLB_a(k,ir) = CRB(1,1)+CRB(2,2);   % angle CRLB
            CRLB_g(k,ir) = CRB(3,3);            % g CRLB
        else
            CRLB_a(k,ir) = CRB(1,1);            % angle CRLB
            CRLB_g(k,ir) = CRB(2,2);            % g CRLB
        end;
        if (N==3)
            E = [-srcLoc(2), -cos(theta)*srcLoc(3), -r(1)*srcLoc(1);
                srcLoc(1), -sin(theta)*srcLoc(3), -r(1)*srcLoc(2);
                0, r(1)*cos(phi), -r(1)*srcLoc(3)];
        else
            E = [ -srcLoc(2), -r(1)*srcLoc(1);
                srcLoc(1), -r(1)*srcLoc(2)];
        end;
        
        CRBp = E*CRB*E';
        CRLB_p(k,ir) = trace(CRBp);         % Cartesian position CRLB
        CRLB_r(k,ir) = srcLoc'*CRBp*srcLoc/(r(1)^2); % source range CRLB
        
        % Theoretical covariance
        Var1 = Cov_SUM_MPR( senPos, srcLoc, Q );
        var_T1(k,ir) = Var1(1);
        var_P1(k,ir) = Var1(2);
        var_g1(k,ir) = Var1(3);
        
        Var2 = Cov_GTRS_MPR( senPos, srcLoc, Q );
        var_T2(k,ir) = Var2(1);
        var_P2(k,ir) = Var2(2);
        var_g2(k,ir) = Var2(3);
        
        % Theoretical bias
        Bias1 = Bias_SUM_MPR( senPos, srcLoc, Q );
        Bias_a1(k,ir) = norm(Bias1(1:end-1));
        Bias_g1(k,ir) = abs(Bias1(end));
        Bias_xy1(k,ir) = norm(E*Bias1,2);
        Bias_Th1(k,ir) = abs(Bias1(1));
        
        Bias2 = Bias_GTRS_MPR( senPos, srcLoc, Q );
        Bias_a2(k,ir) = norm(Bias2(1:end-1));
        Bias_g2(k,ir) = abs(Bias2(end));
        Bias_xy2(k,ir) = norm(E*Bias2,2);
        
        % -- Obtaining source location estimate in MPR --
        nsePwr = 10^(sigmaSquareDB(k)/10);
        rng('default');
        
        [bia_p1,bia_p2,bia_p3,bia_p4,bia_p5,bia_p6] = deal(zeros(N,NumEnsembles));
        for i = 1:NumEnsembles,
            
            % measured TDOAs
            tmp=randn(M,1);
            rdNse = sqrt(nsePwr) * ((tmp(2:M)-tmp(1))/sqrt(2)-PP);
            rd_m = rd + rdNse;
            
            %% SUM-MPR Method
            if (N==3)
                [Th, Ph, g, pos] = TDOA_SUM_MPR( senPos, rd_m, Q );
            else
                [Th, g, pos] = TDOA_SUM_MPR( senPos, rd_m, Q );
                Ph=0;
            end;
            ep1(k,i) = sum(abs(pos-srcLoc).^2);
            eTh1(k,i) = abs(theta - Th)^2;
            ePh1(k,i) = abs(phi - Ph)^2;
            eg1(k,i) = (1/r(1) - g)^2;
            bia_p1(:,i) = pos-srcLoc;
            er1(k,i) = (r(1) - 1/g)^2;
            uTh1(k,i) = Th;
            uPh1(k,i) = Ph;
            ug1(k,i) = g;
            
            %% GTRS-MPR Method
            if (N==3)
                [Th, Ph, g, pos] = TDOA_GTRS_MPR( senPos, rd_m, Q );
            else
                [Th, g, pos] = TDOA_GTRS_MPR( senPos, rd_m, Q );
                Ph=0;
            end;
            ep6(k,i) = sum((pos-srcLoc).^2);
            eTh6(k,i) = (theta - Th)^2;
            ePh6(k,i) = (phi - Ph)^2;
            eg6(k,i) = (1/r(1) - g)^2;
            bia_p6(:,i) = pos-srcLoc;
            er6(k,i) = (r(1) - 1/g)^2;
            uTh6(k,i) = Th;
            uPh6(k,i) = Ph;
            ug6(k,i) = g;
            
            %             % Chan-Ho Method, TSP 1994 (need to download TDOALoc)
            %             pos_ch = TDOALoc(senPos,rd_m,Q);
            %             ep2(k,i) = sum((pos_ch-srcLoc).^2);
            % 			Th = atan2(pos_ch(2),pos_ch(1));
            %             if (N==3)
            %     			Ph = atan2(pos_ch(3),norm(pos_ch(1:2),'fro'));
            %             else
            %                 Ph = 0;
            %             end;
            % 			g = 1/norm(pos_ch-senPos(:,1),'fro');
            %             eTh2(k,i) = (theta-Th).^2;
            %             ePh2(k,i) = (phi-Ph).^2;
            %             eg2(k,i) = (1/r(1) - g).^2;
            %             bia_p2(:,i) = pos_ch-srcLoc;
            %             er2(k,i) = (r(1) - 1/g)^2;
            %             uTh2(k,i) = Th;
            %             uPh2(k,i) = Ph;
            %             ug2(k,i) = g;
            %
            %             % MLE-MPR (need to download TDOA_CVXMPR_3D, TDOA_MLEMPR_3D)
            %             upCVX = TDOA_CVXMPR_3D(senPos,rd_m,Q);           % CVX-MPR solution
            %             upMPR = TDOA_MLEMPR_3D(senPos,rd_m,Q,upCVX);       % MLE-MPR solution
            %             if (N==3)
            %                 pos_MPR = [cos(upMPR(1))*cos(upMPR(2));sin(upMPR(1))*cos(upMPR(2));sin(upMPR(2))]/upMPR(3);
            %             else
            %                 pos_MPR = [cos(upMPR(1));sin(upMPR(1))]/upMPR(2);
            %                 upMPR=[upMPR(1);0;upMPR(2)];
            %             end;
            %             ep5(k,i) = sum((pos_MPR-srcLoc).^2);
            %             eg5(k,i) = (1/r(1) - upMPR(end)).^2;
            %             bia_p5(:,i) = pos_MPR-srcLoc;
            %             er5(k,i) = (r(1) - 1/upMPR(end))^2;
            %             eTh5(k,i) = (theta-upMPR(1)).^2;
            %             ePh5(k,i) = (phi-upMPR(2)).^2;
            %             uTh5(k,i) = upMPR(1);
            %             uPh5(k,i) = upMPR(2);
            %             ug5(k,i) = upMPR(3);
            
        end
        avBia_p1(k,ir) = norm(mean(bia_p1,2),2);
        avBia_p2(k,ir) = norm(mean(bia_p2,2),2);
        avBia_p5(k,ir) = norm(mean(bia_p5,2),2);
        avBia_p6(k,ir) = norm(mean(bia_p6,2),2);
    end
    
    % -- calculate MSE --
    % MSE of angle
    mse_a1(:,ir) = mean(eTh1+ePh1,2);
    mse_a2(:,ir) = mean(eTh2+ePh2,2);
    mse_a5(:,ir) = mean(eTh5+ePh5,2);
    mse_a6(:,ir) = mean(eTh6+ePh6,2);
    
    % MSE of g
    mse_g1(:,ir) = mean(eg1,2);
    mse_g2(:,ir) = mean(eg2,2);
    mse_g5(:,ir) = mean(eg5,2);
    mse_g6(:,ir) = mean(eg6,2);
    
    % MSE of position
    mse_p1(:,ir) = mean(ep1,2);
    mse_p2(:,ir) = mean(ep2,2);
    mse_p5(:,ir) = mean(ep5,2);
    mse_p6(:,ir) = mean(ep6,2);
    
    % MSE of r
    mse_r1(:,ir) = mean(er1,2);
    mse_r2(:,ir) = mean(er2,2);
    mse_r5(:,ir) = mean(er5,2);
    mse_r6(:,ir) = mean(er6,2);
    
    % Bias of angle
    avBia_a1(:,ir) = sqrt((mean(abs(uTh1),2)-theta).^2+(mean(abs(uPh1),2)-phi).^2);%mean(sqrt(eTh1+ePh1),2);
    avBia_a2(:,ir) = sqrt((mean(abs(uTh2),2)-theta).^2+(mean(abs(uPh2),2)-phi).^2);
    avBia_a5(:,ir) = sqrt((mean(abs(uTh5),2)-theta).^2+(mean(abs(uPh5),2)-phi).^2);
    avBia_a6(:,ir) = sqrt((mean(abs(uTh6),2)-theta).^2+(mean(abs(uPh6),2)-phi).^2);
    
    % Bias of g
    avBia_g1(:,ir) = abs(mean(ug1,2)-1/r(1));
    avBia_g2(:,ir) = abs(mean(ug2,2)-1/r(1));
    avBia_g5(:,ir) = abs(mean(ug5,2)-1/r(1));
    avBia_g6(:,ir) = abs(mean(ug6,2)-1/r(1));
end

% ----- plot results -----
clear eTh1 eTh2 eTh3 eTh4 eTh5 eTh6  ePh1 ePh2 ePh3 ePh4 ePh5 ePh6  eg1 eg2 eg3 eg4 eg5 eg6  ...
    ep1 ep2 ep3 ep4 ep5 ep6  bia_p1 bia_p2 bia_p3 bia_p4 bia_p5 bia_p6  er1 er2 er3 er4 er5 er6  ...
    uTh1 uTh2 uTh3 uTh4 uTh5 uTh6  uPh1 uPh2 uPh3 uPh4 uPh5 uPh6  ug1 ug2 ug3 ug4 ug5 ug6;
if length(sigmaSquareDB) ~= 1
    
    indR = 1;
    
    % MSE of angle estimate
    figure;
    plot(sigmaSquareDB, 10*log10(mse_a1(:,indR)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(sigmaSquareDB, 10*log10(mse_a6(:,indR)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(sigmaSquareDB, 10*log10(mse_a2(:,indR)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(sigmaSquareDB, 10*log10(mse_a5(:,indR)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(sigmaSquareDB, 10*log10(CRLB_a(:,indR)), '--', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
    
    xlabel('10log(\sigma_n^2)', 'FontSize', 13);
    ylabel('10log(MSE(\theta,\phi)(rad^2))', 'FontSize', 13);
    lgd11 = legend('Show');
    set(lgd11, 'FontSize',11, 'Location', 'Northwest');
    ylim([-110 20]);
    
    % MSE of g estimate
    figure;
    plot(sigmaSquareDB, 10*log10(mse_g1(:,indR)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(sigmaSquareDB, 10*log10(mse_g6(:,indR)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(sigmaSquareDB, 10*log10(mse_g2(:,indR)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(sigmaSquareDB, 10*log10(mse_g5(:,indR)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(sigmaSquareDB, 10*log10(CRLB_g(:,indR)), '--', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
    
    xlabel('10log(\sigma_n^2)', 'FontSize', 13);
    ylabel('10log(MSE(g)(1/m^2))', 'FontSize', 13);
    lgd2 = legend('Show');
    set(lgd2, 'FontSize',11, 'Location', 'Northwest');
    ylim([-170 -10]);set(gca,'YTick',-170:20:-10);
    
    % MSE of position estimate in Cartesian
    figure;
    plot(sigmaSquareDB, 10*log10(mse_p1(:,indR)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(sigmaSquareDB, 10*log10(mse_p6(:,indR)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(sigmaSquareDB, 10*log10(mse_p2(:,indR)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(sigmaSquareDB, 10*log10(mse_p5(:,indR)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(sigmaSquareDB, 10*log10(CRLB_p(:,indR)), '--', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
    
    xlabel('10log(\sigma_n^2)', 'FontSize', 13);
    ylabel('10log(MSE(u)(m^2))', 'FontSize', 13);
    lgd3 = legend('Show');
    set(lgd3, 'FontSize',11, 'Location', 'Northwest');
    
    % MSE of r estimate
    figure;
    plot(sigmaSquareDB, 10*log10(mse_r1(:,indR)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(sigmaSquareDB, 10*log10(mse_r6(:,indR)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(sigmaSquareDB, 10*log10(mse_r2(:,indR)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(sigmaSquareDB, 10*log10(mse_r5(:,indR)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(sigmaSquareDB, 10*log10(CRLB_r(:,indR)), '--', 'LineWidth', 1.5, 'DisplayName', 'CRLB');%hold on;grid on;
    
    xlabel('10log(\sigma_n^2)', 'FontSize', 13);
    ylabel('10log(MSE(r)(m^{2}))', 'FontSize', 13);
    h1 = legend('Show');
    set(h1, 'FontSize',11, 'Location', 'Northwest');
    
    % Bias of angle estimate
    figure;
    plot(sigmaSquareDB, 20*log10(avBia_a1(:,indR)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(sigmaSquareDB, 20*log10(avBia_a6(:,indR)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(sigmaSquareDB, 20*log10(avBia_a2(:,indR)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(sigmaSquareDB, 20*log10(avBia_a5(:,indR)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(sigmaSquareDB, 20*log10(Bias_a1(:,indR)), '-', 'LineWidth', 1.5, 'DisplayName', 'Thy-SUM','Color',clor(1,:));
    plot(sigmaSquareDB, 20*log10(Bias_a2(:,indR)), '--', 'LineWidth', 1.5, 'DisplayName', 'Thy-GTRS','Color',clor(2,:));
    
    xlabel('10log(\sigma_n^2)', 'FontSize', 13);
    ylabel('20log(Bias(\theta,\phi)(rad))', 'FontSize', 13);
    h3 = legend('Show');
    set(h3, 'FontSize',11, 'Location', 'Northwest');
    xlim([min(sigmaSquareDB),max(sigmaSquareDB)]); ylim([-220 10]);
    
    % Bias of g estimate
    figure;
    plot(sigmaSquareDB, 20*log10(avBia_g1(:,indR)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(sigmaSquareDB, 20*log10(avBia_g6(:,indR)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(sigmaSquareDB, 20*log10(avBia_g2(:,indR)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(sigmaSquareDB, 20*log10(avBia_g5(:,indR)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(sigmaSquareDB, 20*log10(Bias_g1(:,indR)), '-', 'LineWidth', 1.5, 'DisplayName', 'Thy-SUM','Color',clor(1,:));
    plot(sigmaSquareDB, 20*log10(Bias_g2(:,indR)), '--', 'LineWidth', 1.5, 'DisplayName', 'Thy-GTRS','Color',clor(2,:));
    
    xlabel('10log(\sigma_n^2)', 'FontSize', 13);
    ylabel('20log(Bias(g)(m^{-1}))', 'FontSize', 13);
    h3 = legend('Show');
    set(h3, 'FontSize',11, 'Location', 'Northwest');
    xlim([min(sigmaSquareDB),max(sigmaSquareDB)]);ylim([-270 0]);
    
    % Bias of position estimate in Cartesian
    figure;
    plot(sigmaSquareDB, 20*log10(avBia_p1(:,indR)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(sigmaSquareDB, 20*log10(avBia_p1(:,indR)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(sigmaSquareDB, 20*log10(avBia_p2(:,indR)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(sigmaSquareDB, 20*log10(avBia_p5(:,indR)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(sigmaSquareDB, 20*log10(Bias_xy1(:,indR)), '-', 'LineWidth', 1.5, 'DisplayName', 'Thy-SUM','Color',clor(1,:));
    plot(sigmaSquareDB, 20*log10(Bias_xy2(:,indR)), '--', 'LineWidth', 1.5, 'DisplayName', 'Thy-GTRS','Color',clor(2,:));
    
    xlabel('10log(\sigma_n^2)', 'FontSize', 13);
    ylabel('20log(Bias(u)(m))', 'FontSize', 13);
    h3 = legend('Show');
    set(h3, 'FontSize',11, 'Location', 'Southeast');
    xlim([min(sigmaSquareDB),max(sigmaSquareDB)]);ylim([-100,120]);
end

if length(range) ~= 1
    
    % Plot Figure (sigma^2 = 0dB)
    indS = 1;       %find(sigmaSquareDB==0);
    
    % MSE of estimated angle
    figure;
    plot(range/1e3, 10*log10(mse_a1(indS,:)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(range/1e3, 10*log10(mse_a1(indS,:)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(range/1e3, 10*log10(mse_a2(indS,:)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(range/1e3, 10*log10(mse_a5(indS,:)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(range/1e3, 10*log10(CRLB_a(indS,:)), '--', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
    ylim([-70 0]);
    xlabel('Range(km)', 'FontSize', 13);
    ylabel('10log(MSE(\theta,\phi)(rad^2))', 'FontSize', 13);
    lgd11 = legend('Show');
    set(lgd11, 'FontSize',11, 'Location', 'Northwest');
    xlim([min(range),max(range)]/1e3);%ylim([-64.5,-63]);
    
    % MSE of g estimate
    figure;
    plot(range/1e3, 10*log10(mse_g1(indS,:)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(range/1e3, 10*log10(mse_g6(indS,:)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR')
%     plot(range/1e3, 10*log10(mse_g2(indS,:)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(range/1e3, 10*log10(mse_g5(indS,:)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(range/1e3, 10*log10(CRLB_g(indS,:)), '--', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
    ylim([-140 -40]);
    xlabel('Range(km)', 'FontSize', 13);
    ylabel('10log(MSE(g)(1/m^2))', 'FontSize', 13);
    lgd13 = legend('Show');
    set(lgd13, 'FontSize',11, 'Location', 'Northwest');
    xlim([min(range),max(range)]/1e3);%ylim([-120,-118]);
    
    
    % Bias of angle estimate
    figure;
    plot(range/1e3, 20*log10(avBia_a1(indS,:)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(range/1e3, 20*log10(avBia_a6(indS,:)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(range/1e3, 20*log10(avBia_a2(indS,:)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(range/1e3, 20*log10(avBia_a5(indS,:)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(range/1e3, 20*log10(Bias_a1(indS,:)), '-', 'LineWidth', 1.5, 'DisplayName', 'Thy-SUM','Color',clor(1,:));
    plot(range/1e3, 20*log10(Bias_a2(indS,:)), '--', 'LineWidth', 1.5, 'DisplayName', 'Thy-GTRS','Color',clor(2,:));
    ylim([-140 -20]);
    xlabel('Range(km)', 'FontSize', 13);
    ylabel('20log(Bias(\theta,\phi)(rad))', 'FontSize', 13);
    h3 = legend('Show');
    set(h3, 'FontSize',11, 'Location', 'Southeast');
    xlim([min(range/1e3),max(range/1e3)]);
    
    % Bias of g estimate
    figure;
    plot(range/1e3, 20*log10(avBia_g1(indS,:)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');hold on;grid on;
    plot(range/1e3, 20*log10(avBia_g6(indS,:)), 'v', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
%     plot(range/1e3, 20*log10(avBia_g2(indS,:)), '*', 'LineWidth', 1.5, 'DisplayName', 'CFS');
%     plot(range/1e3, 20*log10(avBia_g5(indS,:)), 'x', 'LineWidth', 1.5, 'DisplayName', 'MLE-MPR');
    plot(range/1e3, 20*log10(Bias_g1(indS,:)), '-', 'LineWidth', 1.5, 'DisplayName', 'Thy-SUM','Color',clor(1,:));
    plot(range/1e3, 20*log10(Bias_g2(indS,:)), '--', 'LineWidth', 1.5, 'DisplayName', 'Thy-GTRS','Color',clor(2,:));
    ylim([-200 -40]);
    xlabel('Range(km)', 'FontSize', 13);
    ylabel('20log(Bias(g)(m^{-1}))', 'FontSize', 13);
    h3 = legend('Show');
    set(h3, 'FontSize',11, 'Location', 'Southeast');
    xlim([min(range/1e3),max(range/1e3)]);%ylim([-270 0]);
    
end


