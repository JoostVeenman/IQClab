% -------------------------------------------------------------------------
%
% IQClab:      Version 3.04
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NoDerivatives 4.0
%              International (CC BY-ND 4.0)) license:  
%              https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        17-01-2020
% 
% -------------------------------------------------------------------------
% Demo_006:    Robustness analysis with sector-bounded and slope restricted
%              scalar nonlinearities
%
%              This demo file performs an IQC-robustness analysis for an
%              uncertain plant that is affected by sector-bounded and slope
%              restricted scalar nonlinearities. Here it is possible to
%              specify various options: 
%
%                1.) Uncertainty block  
%                      a.) 3x repeated sector-bounded scalar nonlinearity
%                      b.) 3 different sector-bounded scalar nonlinearities
%                      c.) 3x repeated sector-bounded and slope restricted
%                          scalar nonlinearity
%                      d.) 3 different sector-bounded and slope-restricted
%                          scalar nonlinearities 
%
% -------------------------------------------------------------------------
clc;close all;clear ga

% Specify inputs
clc
disp('------------------------------------------------------------------');
disp('Demo_006: Robustness analysis with sector-bounded and slope');
disp('          restricted scalar nonlinearities');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file performs an IQC robustness analysis for an');
disp('uncertain plant that is affected by sector-bounded and slope');
disp('restricted scalar nonlinearities. Here it is possible to specify');
disp('the following options:');
disp(' ');
disp('  1.) Uncertainty block');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Uncertainty block - Please choose from the following options:');
disp('  1: 3 different sector-bounded nonlinearities, each repeated once');
disp('  2: 3x repeated sector-bounded nonlinearity');
disp('  3: 3 different sector-bounded and slope restricted');
disp('     nonlinearities, each repeated once')
disp('  4: 3x repeated sector-bounded and slope restricted nonlinearity');
disp(' ');
nu = input('Choose 1-4: ');
disp(' ');
disp('------------------------------------------------------------------');

% Define plant
A                       = [-0.4,-1;1,0];
Bp                      = [-0.2,-1,-0.25;0,0,0];
Bw                      = [1;0];
Cq                      = [1,0;0,1;0,0];
Cz                      = [0,-0.2];
Dqp                     = [0,0,0;0,0,0;0,1,0];
Dqw                     = [0;0;0];
Dzp                     = [-0.1,0,0];
Dzw                     = 1;
N                       = ss(A,[Bp,Bw],[Cq;Cz],[Dqp,Dqw;Dzp,Dzw]);


switch nu
    case 2
        % Define uncertainty range
        alpha           = linspace(0,0.7,25)';alpha(1) = 0.00001;
        l               = 1;
    case 1
        % Define uncertainty range
        alpha           = linspace(0,0.7,25)';alpha(1) = 0.00001;
        l               = 1;
    case 4
        % Define uncertainty range
        alpha1          = linspace(0,0.9,15)';alpha1(1) = 0.00001;
        alpha2          = linspace(0,2.25,15)';alpha2(1) = 0.00001;
        alpha3          = linspace(0,2.25,15)';alpha3(1) = 0.00001;
        alpha4          = linspace(0,5.5,15)';alpha4(1) = 0.00001;
        alpha           = [alpha1,alpha2,alpha3,alpha4];
        l               = 1:1:4;
    case 3
        % Define uncertainty range
        alpha1          = linspace(0,0.7,15)';alpha1(1) = 0.00001;
        alpha2          = linspace(0,1.5,15)';alpha2(1) = 0.00001;
        alpha3          = linspace(0,1.5,15)';alpha3(1) = 0.00001;
        alpha4          = linspace(0,3.7,15)';alpha4(1) = 0.00001;
        alpha           = [alpha1,alpha2,alpha3,alpha4];
        l               = 1:1:4;
end

for j = 1:length(l)
    for i = 1:length(alpha)
        switch nu
            case 2
                delta   = iqcdelta('delta','InputChannel',1:3,'OutputChannel',1:3,'LinNonlin','NL','SectorBounds',[0,alpha(i,j)]);
                udelta  = iqcassign(delta,'usbsr');
            case 1
                delta1  = iqcdelta('delta1','InputChannel',1,'OutputChannel',1,'LinNonlin','NL','SectorBounds',[0,alpha(i,j)]);
                delta2  = iqcdelta('delta2','InputChannel',2,'OutputChannel',2,'LinNonlin','NL','SectorBounds',[0,alpha(i,j)]);
                delta3  = iqcdelta('delta3','InputChannel',3,'OutputChannel',3,'LinNonlin','NL','SectorBounds',[0,alpha(i,j)]);
                udelta1 = iqcassign(delta1,'usbsr');
                udelta2 = iqcassign(delta2,'usbsr');
                udelta3 = iqcassign(delta3,'usbsr');
            case 4
                delta   = iqcdelta('delta','InputChannel',1:3,'OutputChannel',1:3,'LinNonlin','NL','Odd','yes','SectorBounds',[0,alpha(i,j)],'SlopeBounds',[0,alpha(i,j)]);
                udelta  = iqcassign(delta,'usbsr','Length',l(j),'PoleLocation',-1);
            case 3
                delta1  = iqcdelta('delta1','InputChannel',1,'OutputChannel',1,'LinNonlin','NL','Odd','yes','SectorBounds',[0,alpha(i,j)],'SlopeBounds',[0,alpha(i,j)]);
                delta2  = iqcdelta('delta2','InputChannel',2,'OutputChannel',2,'LinNonlin','NL','Odd','yes','SectorBounds',[0,alpha(i,j)],'SlopeBounds',[0,alpha(i,j)]);
                delta3  = iqcdelta('delta3','InputChannel',3,'OutputChannel',3,'LinNonlin','NL','Odd','yes','SectorBounds',[0,alpha(i,j)],'SlopeBounds',[0,alpha(i,j)]);
                udelta1 = iqcassign(delta1,'usbsr','Length',l(j),'PoleLocation',-1);
                udelta2 = iqcassign(delta2,'usbsr','Length',l(j),'PoleLocation',-1);
                udelta3 = iqcassign(delta3,'usbsr','Length',l(j),'PoleLocation',-1);
        end
        perf            = iqcdelta('perf','ChannelClass','P','InputChannel',4,'OutputChannel',4,'PerfMetric','L2');
        switch nu
            case {1,3}
                Delta   = {udelta1,udelta2,udelta3,perf};
            case {2,4}
                Delta   = {udelta,perf};
        end
        prob            = iqcanalysis(N,Delta);
        ga(i,j)         = prob.gamma;
        disp(prob.gamma);
    end
end
% Plot results
lt                      = {'*r--','ob:','<g-.','xk-'};
for j = 1:length(l)
    se                  = ga(:,j) > 0;
    figure(1)
    plot(alpha(se,j),ga(se,j),lt{j});hold on
end
switch nu
    case {1,2}
        axis([0,alpha(end,j),1,3]);
        xlabel('\alpha, \delta\in sector(0,\alpha)');
    case {3,4}
        axis([0,alpha(end,j),1,3]);
        xlabel('\alpha, \delta\in sector(0,\alpha)\cap slope(0,\alpha)');
        legend('iqc, l=1','iqc, l=2','iqc, l=3','iqc, l=4','Location','NE');
end
ylabel('\gamma');
grid on
title('Worst-case induced L_2-gain');
set(gcf,'Position',[300,400,700,350]);
fCutFig(1,1);