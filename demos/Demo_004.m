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
% Demo_004:    Robustness analysis with arbitrarily fast varying parametric
%              uncertainties
%
%              This demo file performs an IQC robustness analysis for an
%              uncertain plant that is affected by arbitrarily fast varying
%              parametric uncertainties and using different relaxation
%              schemes. Here it is possible to specify various options: 
%
%                1.) Uncertainty block (two parametric uncertainties that
%                    are repeated once, or four LTV parametric
%                    uncertainties packed in a full 2x2-block
%                2.) Performance metric (L2, H2, robust stability only)
%
% -------------------------------------------------------------------------
clc;close all;clear ga

% Specify inputs
clc
disp('------------------------------------------------------------------');
disp('Demo_004: Robustness analysis with arbitrarily fast varying');
disp('          parametric uncertainties');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file performs an IQC robustness analysis for an');
disp('uncertain plant that is affected by arbitrarily fast varying.');
disp('parametric uncertainties and using different relaxation schemes.');
disp('Here it is possible to specify the following options:');
disp(' ');
disp('  1.) Uncertainty block');
disp('  2.) Performance metric');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Uncertainty block - Please choose from the following options:');
disp('  1: Two LTV parametric uncertainties, each repeated once');
disp('     Delta = blkdiag(delta1,delta2)');
disp('  2: Four LTV parametric uncertainties, 2x2 full-block');
disp('     Delta = H1*delta1 + ... + H4*delta4'); 
disp(' ');
nu = input('Choose 1 or 2: ');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Performance metric - Please choose from the following options:');
disp('  1: Worst-case induced L2-gain');
disp('  2: Worst-case H2-norm');
disp('  3: Robust stability test');
disp(' ');
p = input('Choose 1, 2 or 3: ');
disp(' ');
disp('------------------------------------------------------------------');


% Define plant
A                               = [-2,-3;1,1];
Bp                              = [1,0;0,0];
Bw                              = [1;0];
Cq                              = [1,0;0,0];
Cz                              = [1,0];
Dqp                             = [1,-2;1,-1];
switch p
    case {1,3}
        Dqw                     = [0;1];
    case 2
        Dqw                     = [0;0];
end
Dzp                             = [0,1];
Dzw                             = 0;
N                               = ss(A,[Bp,Bw],[Cq;Cz],[Dqp,Dqw;Dzp,Dzw]);

% Define uncertainty range
switch nu
    case 1
        alpha                   = linspace(0,0.4,50);alpha(1) = 0.0001;
        rt                      = {'DG','CH','PC','ZP'};
    case 2
        alpha                   = linspace(0,0.2,50);alpha(1) = 0.0001;
        rt                      = {'CH','ZP'};
end
for j = 1:length(rt)
    for i = 1:length(alpha)
        switch nu
            case 1
                H1{1}           = [1,0;0,0];
                H1{2}           = [0,0;0,1];
                La              = polydec(pvec('box',alpha(i)*[-1,1;-1,1]))';

                delta           = iqcdelta('delta','Structure','FB','UncertaintyMap',H1,'Polytope',La,'InputChannel',1:2,'OutputChannel',1:2,'TimeInvTimeVar','TV');
                udelta          = iqcassign(delta,'ultv','RelaxationType',rt{j});
            case 2
                H2{1}           = [1,0;0,1];
                H2{2}           = [1,1;0,0];
                H2{3}           = [0,1;1,0];
                H2{4}           = [1,0;1,0];
                La              = polydec(pvec('box',alpha(i)*[-1,1;-1,1;-1,1;-1,1]))';
                delta           = iqcdelta('delta','Structure','FB','UncertaintyMap',H2,'Polytope',La,'InputChannel',1:2,'OutputChannel',1:2,'TimeInvTimeVar','TV');
                udelta          = iqcassign(delta,'ultv','RelaxationType',rt{j});
        end
        switch p
                case 1
                    perf        = iqcdelta('perf','ChannelClass','P','InputChannel',3,'OutputChannel',3,'PerfMetric','L2');
                    Delta       = {udelta,perf};
                case 2
                    perf        = iqcdelta('perf','ChannelClass','P','InputChannel',3,'OutputChannel',3,'PerfMetric','H2');
                    Delta       = {udelta,perf};
                case 3
                    Delta       = udelta;
        end
        prob                    = iqcanalysis(N,Delta);
        switch p
            case {1,2}
                ga(i,j)         = prob.gamma;
                disp(prob.gamma);
            case 3
                if strcmp(prob.gamma,'feasible')
                    ga(i,j)         = 1;
                else
                    ga(i,j) = -1;
                end
                disp(prob.gamma);
        end
    end
end
% Plot results
lt = {'+b-','*r--','oc:','<g-.','xk-'};
switch p
    case {1,2}
        figure(1)
        for j = 1:length(rt)
            se = ga(:,j) > 0;
            figure(1)
            plot(alpha(se)',ga(se,j),lt{j+1});hold on
        end
        switch nu
            case 1
                xlabel('\alpha (\delta_i\in[-\alpha,\alpha], i = 1,2)');
                legend('RelaxType, DG','RelaxType, CH','RelaxType, PC','RelaxType, ZP');
            case 2
                xlabel('\alpha (\delta_i\in[-\alpha,\alpha], i = 1,...,4)');
                legend('RelaxType, CH','RelaxType, ZP');
        end
        axis([0,alpha(end),0.5,10]);
        ylabel('\gamma');
        switch p
            case 1
                title('Worst-case induced L_2-gain');
            case 2
                title('Worst-case H_2-norm');
        end
    case 3
        figure(1)
        for j = 1:length(rt)
            figure(1)
            plot(alpha',ga(:,j),lt{j+1});hold on
        end
        switch nu
            case 1
                xlabel('\alpha (\delta_i\in[-\alpha,\alpha], i = 1,2)');
                legend('RelaxType, DG','RelaxType, CH','RelaxType, PC','RelaxType, ZP','Location','SW');
            case 2
                xlabel('\alpha (\delta_i\in[-\alpha,\alpha], i = 1,...,4)');
                legend('RelaxType, CH','RelaxType, ZP','Location','SW');
        end
        axis([0,alpha(end),0.5,10]);
        ylabel('1 = rob. stable, -1 = LMIs infeasible');
        axis([0,alpha(end),-1.1,1.1]);
        title('Robust stabilty');
end
grid on
set(gcf,'Position',[300,400,700,350]);
fCutFig(1,1);