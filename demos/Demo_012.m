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
% Date:        07-04-2020
% 
% -------------------------------------------------------------------------
% Demo_012:    Control mode switching
%
%              This demo demonstrates the potential effectiveness of the
%              Youla-based control switching implementation scheme.
% -------------------------------------------------------------------------
close all;clear all;clc;set(0,'DefaultAxesFontSize',14);

disp('------------------------------------------------------------------');
disp('Demo_012: Control mode switching');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file compares the performance of different switching');
disp('schemes. In particular, the demo demonstrates the potential');
disp('benefit of implementing the so called Youla based switching scheme');
disp(' ');
disp('Two different switching implementations are compared:');
disp('  1.) The classical implementation: K = (1-alpha)*K1 + alpha K2');
disp('  2.) The Youla-based implementation: Q = (1-alpha)*Q1 + alpha Q2');
disp('Here K1 and K2 are two controller and alpha\in[0,1] is the');
disp('switching parameter, while Q1 and Q2 are Youla parameters obtained')
disp('from the controllers K1 and K2 respectively.');
disp(' ');
disp('It is possible to consider two different switching schemes:');
disp('  1. Instantaneous switching: alpha(t) is instantaneously switched');
disp('     from 0 to 1.')
disp('  2.) Smooth switching: alpha(t) is slowly increased from 0 to 1.');
disp(' ');
disp('------------------------------------------------------------------');
disp(' ');
disp('Please choose from the following options:');
disp('  0: Smooth switch');
disp('  1: Instanteneous switch');
disp(' ');
nu = input('Choose 0 or 1: ');
disp(' ');
disp('------------------------------------------------------------------');

mode            = nu; % 0 = smooth switch, 1 = instanteneous switch
N               = 20;
bw              = ureal('bw',5,'Percentage',10);
Gnom            = tf(1,[1/bw,1]);
W               = makeweight(.05,9,10);
Delta           = ultidyn('Delta',[1 1]);
G               = Gnom*(1+W*Delta);

systemnames     = 'G';
inputvar        = '[r;n;u]';
outputvar       = '[r-G;u;G;r-G-n]';
input_to_G      = '[u]';
cleanupsysic    = 'yes';
uolic           = sysic;
Usamp           = usample(uolic,N);
Usamp(:,:,1,1)  = uolic.NominalValue;
Usamp(:,:,N,1) = uolic.NominalValue;

xi              = 0.707;
wn1             = 3;
wn2             = 7.5;

Kp1             = 2*xi*wn1/5 - 1;
Ki1             = (wn1^2)/5;
K1              = ss(tf([Kp1,Ki1],[1 0]));

Kp2             = 2*xi*wn2/5 - 1;
Ki2             = (wn2^2)/5;
K2              = ss(tf([Kp2,Ki2],[1 0]));

for i = 1:N
    olic        = Usamp(:,:,i,1);
    G22         = olic(4,3);
    [L,Q]       = fYoulaSwitch(G22,K1,K2);
    Q1          = Q{1};
    Q2          = Q{2};
    [nQo,nQi]   = size(Q1.d);
    [nKo,nKi]   = size(K1.d);

    sim('YoulaSwitch');

    figure(1)
    subplot(1,2,1)
    if i == 1
        plot(disturbances.Time,disturbances.Data(:,1),'k--','LineWidth',1);hold on
        plot(plant_output.Time,plant_output.Data(:,1),'k-','LineWidth',2);
    elseif i == N
        plot(plant_output.Time,plant_output.Data(:,1),'k-','LineWidth',2);
    else
        plot(plant_output.Time,plant_output.Data(:,1),'r:','LineWidth',1);
    end

    subplot(1,2,2)
    if i == 1
        plot(disturbances.Time,disturbances.Data(:,1),'k--','LineWidth',1);hold on
        plot(plant_output.Time,plant_output.Data(:,2),'k-','LineWidth',2);
    elseif i == N
        plot(plant_output.Time,plant_output.Data(:,2),'k-','LineWidth',2);
    else
        plot(plant_output.Time,plant_output.Data(:,2),'r:','LineWidth',1);
    end
end

figure(1)
subplot(1,2,1)
if mode == 0
    axis([0,10,-0.4,1.2]);
elseif mode == 1
    axis([0,10,-1.5,5]);
end
title('Youla-based switch');
grid on
legend('Reference Signal','Nominal response','Uncertain response','Location','SW');

subplot(1,2,2)
if mode == 0
    axis([0,10,-0.4,1.2]);
elseif mode == 1
    axis([0,10,-1.5,5]);
end
title('Classical switch');
grid on
set(gcf,'Position',[23,44,1321,640]);
fCutFig(1,2)