% -------------------------------------------------------------------------
%
% IQClab:      Version 3.4.0
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
% Demo_005:    Robustness analysis with rate-bounded time-varying
%              parametric uncertainties
%
%              This demo file performs an IQC robustness analysis for an
%              uncertain plant that is affected by rate-bounded
%              time-varying parametric uncertainties and using different
%              relaxation schemes. Here it is possible to specify various
%              options:
%
%                1.) Uncertainty block (One rate-bounded LTV parametric
%                    uncertainty which is repeated twice or two
%                    rate-bounded LTV parametric uncertainties which are
%                    repeated once. 
%                2.) Relaxation type ('DG', 'CH', 'PC', 'ZP').
%
% -------------------------------------------------------------------------
clc;close all;clear ga

% Specify inputs
clc
disp('------------------------------------------------------------------');
disp('Demo_005: Robustness analysis with rate-bounded time-varying');
disp('          parametric uncertainties');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file performs an IQC robustness analysis for an');
disp('uncertain plant that is affected by rate-bounded time-varying');
disp('parametric uncertainties and using different relaxation schemes.');
disp('Here it is possible to specify the following options:');
disp(' ');
disp('  1.) Uncertainty block');
disp('  2.) Relaxation type');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Uncertainty block - Please choose from the following options:');
disp('  1: One LTV rate-bounded parameter, repeated twice');
disp('  2: Two LTV rate-bounded parameters, each repeated once');
disp(' ');
nu = input('Choose 1 or 2: ');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Relaxation type - Please choose from the following options:');
disp('  DG: Relaxation type DG scalings');
disp('  CH: Relaxation type convex hull');
disp('  PC: Relaxation type pavrtial convexity');
disp('  ZP: Relaxation type zeroth order Polya');
disp(' ');
rt = input('Choose DG, CH, PC, ZP: ','s');
disp(' ');
disp('------------------------------------------------------------------');

% Define plant
A                      = [-2,-3;1,1];
Bp                     = [1,0;0,0];
Bw                     = [1;0];
Cq                     = [1,0;0,0];
Cz                     = [1,0];
Dqp                    = [1,-2;1,-1];
Dqw                    = [0;1];
Dzp                    = [0,1];
Dzw                    = 0;
N                      = ss(A,[Bp,Bw],[Cq;Cz],[Dqp,Dqw;Dzp,Dzw]);

% Define uncertainty range
alpha                  = 0.3;
alphadot               = linspace(0,2,15);alphadot(1) = 0.00001;

% Perform mu-analysis
switch nu
    case 1
        de             = ureal('de',0,'Range',alpha*[-1,1])*eye(2);        
    case 2
        de1            = ureal('de1',0,'Range',alpha*[-1,1]);
        de2            = ureal('de2',0,'Range',alpha*[-1,1]);
        de             = blkdiag(de1,de2);
end
Mcl                    = lft(de,N);
[wcg,wcu]              = wcgain(Mcl);
if wcg.UpperBound == Inf
    gam                = -1;
else
    gam                = wcg.UpperBound;
end
disp(gam);

% Perform IQC-analysis
l                      = 1:1:3;
for j = 1:length(l)
    for i = 1:length(alphadot)
        switch nu
            case 1
                delta  = iqcdelta('delta','TimeInvTimeVar','TV','InputChannel',1:2,'OutputChannel',1:2,'Bounds',alpha*[-1,1],'RateBounds',alphadot(i)*[-1,1]);
                udelta = iqcassign(delta,'ultv_rb','Length',l(j),'RelaxationType',rt);
            case 2
                delta1 = iqcdelta('delta1','TimeInvTimeVar','TV','InputChannel',1,'OutputChannel',1,'Bounds',alpha*[-1,1],'RateBounds',alphadot(i)*[-1,1]);
                delta2 = iqcdelta('delta2','TimeInvTimeVar','TV','InputChannel',2,'OutputChannel',2,'Bounds',alpha*[-1,1],'RateBounds',alphadot(i)*[-1,1]);
                delta  = blkdiag('delta',delta1,delta2);
                udelta = iqcassign(delta,'ultv_rb','Length',l(j),'RelaxationType',rt);
        end
        perf           = iqcdelta('perf','ChannelClass','P','InputChannel',3,'OutputChannel',3,'PerfMetric','L2');
        Delta          = {udelta,perf};
        prob           = iqcanalysis(N,Delta);
        ga(i,j)        = prob.gamma;
        disp(prob.gamma);
    end
end
% Plot results
lt = {'+b-','*r--','oc:','xk-.'};
figure(1)
plot([0,alphadot(end)],gam*[1,1],lt{1});hold on
for j = 1:length(l)
    se = ga(:,j) > 0;
    figure(1)
    plot(alphadot(se)',ga(se,j),lt{j+1});hold on
end
switch nu
    case 1
        text = ['alphadot (deltadot\in alphadot[-1,1], delta\in[-',num2str(alpha),',',num2str(alpha),'])'];
        xlabel(text);
    case 2
        text = ['alphadot (deltadoti\in alphadot[-1,1], deltai\in[-',num2str(alpha),',',num2str(alpha),'], i = 1,2)'];
        xlabel(text);
end
ylabel('\gamma');
grid on
legend('\mu','iqc, l=1','iqc, l=2','iqc, l=3','Location','SE');
set(gcf,'Position',[300,400,700,350]);
fCutFig(1,1);