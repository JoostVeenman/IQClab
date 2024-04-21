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
% Demo_002:    Robustness analysis with LTI dynamic uncertainties
%
%              This demo file performs a mu and IQC robustness analysis for
%              an uncertain plant that is affected by LTI dynamic
%              uncertainties. Here it is possible to specify various
%              options:
%
%                1.) Uncertainty block (one LTI dynamic 2x2 full-block
%                    uncertainty or an LTI dynamic 1x1 uncertainty that is
%                    repeated twice.
%                2.) Performance metric (L2, H2, robust stability only)
%
% -------------------------------------------------------------------------
clc;close all;clear ga

% Specify inputs
clc
disp('------------------------------------------------------------------');
disp('Demo_002: Robustness analysis with LTI dynamic uncertainties');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file performs a mu and IQC robustness analysis for an');
disp('uncertain plant that is affected by LTI dynamic uncertainties.');
disp('Here it is possible to specify the following options:');
disp(' ');
disp('  1.) Uncertainty block');
disp('  2.) Performance metric');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Uncertainty block - Please choose from the following options:');
disp('  1: One 2x2 LTI full-block dynamic uncertainty, repeated once');
disp('  2: Two 1x1 LTI scalar dynamic uncertainties, repeated twice');
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
A                           = [-2,-3;1,1];
Bp                          = [1,0;0,0];
Bw                          = [1;0];
Cq                          = [1,0;0,0];
Cz                          = [1,0];
Dqp                         = [1,-2;1,-1];
switch p
    case {1,3}
        Dqw                 = [0;1];
    case 2
        Dqw                 = [0;0];
end
Dzp                         = [0,1];
Dzw                         = 0;
N                           = ss(A,[Bp,Bw],[Cq;Cz],[Dqp,Dqw;Dzp,Dzw]);

% Define uncertainty range
alpha                       = linspace(0,0.5,50);alpha(1) = 0.001;

% Perform mu-analysis
for i = 1:length(alpha)
    switch nu
        case 1
            de              = ultidyn('de',[2,2],'Bound',alpha(i));
        case 2
            de              = ultidyn('de',[1,1],'Bound',alpha(i))*eye(2);
    end
    M                       = N;
    Mcl                     = lft(de,M);
    [wcg,wcu]               = wcgain(Mcl);
    switch p
        case {1,2}
            if wcg.UpperBound == Inf
                gam(i)      = -1;
            else
                gam(i)      = wcg.UpperBound;
            end
            disp(gam(i));
        case 3
            if wcg.UpperBound == Inf
                gam(i)      = -1;
                disp('infeasible');
            else
                gam(i)      = 1;
                disp('feasible');
            end
    end
end

% Perform IQC-analysis
l                           = 1:1:4;
for j = 1:length(l)
    for i = 1:length(alpha)
        switch nu
            case 1
                delta       = iqcdelta('delta','InputChannel',1:2,'OutputChannel',1:2,'StaticDynamic','D','Structure','FB','NormBounds',alpha(i));
                udelta      = iqcassign(delta,'ultid','Length',l(j));
            case 2
                delta       = iqcdelta('delta','InputChannel',[1;2],'OutputChannel',[1;2],'StaticDynamic','D','NormBounds',alpha(i));
                udelta      = iqcassign(delta,'ultid','Length',l(j));
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
        prob                = iqcanalysis(N,Delta);
        switch p
            case {1,2}
                ga(i,j)     = prob.gamma;
                disp(prob.gamma);
             case 3
                if strcmp(prob.gamma,'feasible')
                    ga(i,j) = 1;
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
        se = gam > 0;
        plot(alpha(se),gam(se),lt{1});hold on
        for j = 1:length(l)
            se = ga(:,j) > 0;
            figure(1)
            plot(alpha(se)',ga(se,j),lt{j+1});hold on
        end
        axis([0,alpha(end),0.5,10]);
        xlabel('\alpha (||\delta||_\infty<\alpha)');
        ylabel('\gamma');
        switch p
            case 1
                title('Worst-case induced L_2-gain');
            case 2
                title('Worst-case H_2-norm');
        end
        legend('mu','iqc, l=1','iqc, l=2','iqc, l=3','iqc, l=4');
    case 3
        figure(1)
        plot(alpha,gam,lt{1});hold on
        for j = 1:length(l)
            figure(1)
            plot(alpha',ga(:,j),lt{j+1});hold on
        end
        xlabel('\alpha (||\delta||_\infty<\alpha)');
        ylabel('1 = rob. stable, -1 = LMIs infeasible');
        axis([0,alpha(end),-1.1,1.1]);
        title('Robust stabilty');
        legend('mu','iqc, l=1','iqc, l=2','iqc, l=3','iqc, l=4','Location','SW');
end
grid on
set(gcf,'Position',[300,400,700,350]);
fCutFig(1,1);