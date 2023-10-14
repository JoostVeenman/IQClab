% -------------------------------------------------------------------------
%
% IQClab:      Version 3.03
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NonCommercial-NoDerivatives
%              4.0 International (CC BY-NC-ND 4.0))license: 
%              https://creativecommons.org/licenses/by-nc-nd/4.0/
%              Commercial usage is only permitted with a commercial
%              license. For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        17-01-2020
% 
% -------------------------------------------------------------------------
% Demo_003:    Robustness analysis with LTI delay uncertainties
%
%              This demo file performs a mu and IQC robustness analysis for
%              an uncertain plant that is affected by delay uncertainties.
%              Here it is possible to specify various options:
%
%                1.) xtrIQC (Add extra IQC constraint).
%                2.) Performance metric (L2, robust stability only)
%
% -------------------------------------------------------------------------
clc;close all;clear ga

% Specify inputs
clc
disp('------------------------------------------------------------------');
disp('Demo_003: Robustness analysis with LTI delay uncertainties');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file performs a mu and IQC robustness analysis for an');
disp('uncertain plant that is affected by delay uncertainties. Here it');
disp('is possible to specify the following options:');
disp(' ');
disp('  1.) xtrIQC');
disp('  2.) Performance metric');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('xtrIQC - Please choose from the following options:');
disp('  1: Use the standard IQC for delay uncertainties');
disp('  2: Add extra IQC constraint on top of the standard delay IQC');
disp(' ');
xtrIQC = input('Choose 1 or 2: ');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Performance metric - Please choose from the following options:');
disp('  1: Worst-case induced L2-gain');
disp('  2: Robust stability test');
disp(' ');
p = input('Choose 1 or 2: ');
disp(' ');
disp('------------------------------------------------------------------');

% Define plant
s                           = tf('s');
G                           = ss(4/(s^2+0.1*s+1));
Gd                          = ss(10/(s+0.1));
We                          = 2/3*ss((s+3.674)^2/(s+0.03)^2);
Wu                          = ss((s+10)/(s+1e4));
systemnames                 = 'G Gd';
inputvar                    = '[p;d;r;u]';
outputvar                   = '[Gd;r-p-Gd;u;r-p-Gd]';
input_to_G                  = '[u]';
input_to_Gd                 = '[G+d]';
cleanupsysic                = 'yes';
olic                        = sysic;
wolic                       = ssbal(blkdiag(1,We,Wu,1)*olic);
[K,CL,gamma]                = hinfsyn(wolic(2:end,2:end),1,1);
N                           = minreal(lft(wolic,K));

% Define uncertainty range
alpha                       = linspace(0,0.05,50);alpha(1) = 0.0001;

% Perform mu-analysis
for i = 1:length(alpha)
    de                      = ultidyn('de',[1,1],'Bound',1);
    W                       = fW(alpha(i),1e-6,3);
    M                       = N*blkdiag(W,eye(2));
    Mcl                     = lft(de,M);
    [wcg,wcu]               = wcgain(Mcl);
    switch p
        case 1
            if wcg.UpperBound == inf
                gam(i)      = -1;
            else
                gam(i)      = wcg.UpperBound;
            end
            disp(gam(i));
        case 2
            if wcg.UpperBound == inf
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
        switch xtrIQC
            case 1
                delta       = iqcdelta('delta','InputChannel',1,'OutputChannel',1,'StaticDynamic','D','DelayType',1,'DelayTime',alpha(i));
                udelta      = iqcassign(delta,'udel','Length',l(j),'AddIQC','no');
            case 2
                delta       = iqcdelta('delta','InputChannel',1,'OutputChannel',1,'StaticDynamic','D','DelayType',1,'DelayTime',alpha(i));
                udelta      = iqcassign(delta,'udel','Length',l(j),'AddIQC','yes');
        end
        switch p
            case 1
                perf        = iqcdelta('perf','ChannelClass','P','InputChannel',2:3,'OutputChannel',2:3,'PerfMetric','L2');
                Delta       = {udelta,perf};
            case 2
                Delta       = udelta;
        end
        prob                = iqcanalysis(N,Delta);
        switch p
            case 1
                ga(i,j)     = prob.gamma;
                disp(prob.gamma);
            case 2
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
    case 1
        figure(1)
        se = gam > 0;
        plot(alpha(se),gam(se),lt{1});hold on
        for j = 1:length(l)
            se = ga(:,j) > 0;
            figure(1)
            plot(alpha(se)',ga(se,j),lt{j+1});hold on
        end
        axis([0,alpha(end),0.5,10]);
        xlabel('\alpha (max. delay-time < \alpha)');
        ylabel('\gamma');
        title('Worst-case induced L_2-gain');
        legend('mu','iqc, l=1','iqc, l=2','iqc, l=3','iqc, l=4')
    case 2
        figure(1)
        plot(alpha,gam,lt{1});hold on
        for j = 1:length(l)
            figure(1)
            plot(alpha',ga(:,j),lt{j+1});hold on
        end
        xlabel('\alpha (max. delay-time < \alpha)');
        axis([0,alpha(end),0.5,10]);
        ylabel('1 = rob. stable, -1 = LMIs infeasible');
        axis([0,alpha(end),-1.1,1.1]);
        title('Robust stabilty');
        legend('mu','iqc, l=1','iqc, l=2','iqc, l=3','iqc, l=4','Location','SW');
end
grid on
set(gcf,'Position',[300,400,700,350]);
fCutFig(1,1);