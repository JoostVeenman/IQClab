% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        27-03-2021
% 
% -------------------------------------------------------------------------
% Demo_010:    This demo demonstrates how to design an anti-windup
%              compensator using fAWsyn
%
% -------------------------------------------------------------------------
clc;close all;clc;set(0,'DefaultAxesFontSize',14);

% Specify inputs
disp('------------------------------------------------------------------');
disp('Demo_010: Anti-windup design:');
disp('------------------------------------------------------------------');
disp(' ');

% Define random symmetric matrix of dimension n
s               = tf('s');
G               = ss(10/(100*s+1)*[4,-5;-3,4]);
K               = ss(zeros(2),eye(2),[2,2.5;1.5,2]/100,[2,2.5;1.5,2]);
Ke              = ss(K.a,[K.b,eye(2),zeros(2)],K.c,[K.d,zeros(2),eye(2)]);

systemnames     = 'G Ke';
inputvar        = '[p{2};r{2};u{4}]';
outputvar       = '[Ke;r-G;p]';
input_to_G      = '[Ke-p]';
input_to_Ke     = '[r-G;u]';
cleanupsysic    = 'yes';
Paw             = sysic;

AWopt.FeasbRad  = 1e4;
AWopt.constants = [1e-6,1e-6,1e-6];
AWopt.subopt    = 1.03;
[Kaw,ga]        = fAWsyn(Paw,[2,2,2],[2,2,4],AWopt);

open('AW_comp');
out = sim('AW_comp');
close all
t   = out.A.Time;
r1  = out.A.data(:,1);
r2  = out.A.data(:,2);
u1A = out.A.data(:,7);
u2A = out.A.data(:,8);
u1B = out.B.data(:,7);
u2B = out.B.data(:,8);
u1C = out.C.data(:,7);
u2C = out.C.data(:,8);
y1A = out.A.data(:,9);
y2A = out.A.data(:,10);
y1B = out.B.data(:,9);
y2B = out.B.data(:,10);
y1C = out.C.data(:,9);
y2C = out.C.data(:,10);

figure(1)
title('System response');
subplot(2,1,1)
plot(t,r1,'k-');hold on
plot(t,y1A,'r:','LineWidth',2);
plot(t,y1B,'g--','LineWidth',2);
plot(t,y1C,'b-.','LineWidth',2);
axis([0,200,-6,4]);
grid minor

subplot(2,1,2)
plot(t,r2,'k-');hold on
plot(t,y2A,'r:','LineWidth',2);
plot(t,y2B,'g--','LineWidth',2);
plot(t,y2C,'b-.','LineWidth',2);
axis([0,200,-2,5]);
grid minor
legend('Reference','Nominal response without saturation','Nominal response with saturation','Response with AW compensator');
set(gcf,'Position',[23,44,1321,640]);
fCutFig(2,1)

figure(2)
title('Control commands');
subplot(2,1,1)
plot(t,u1A,'r:','LineWidth',2);hold on
plot(t,u1B,'g--','LineWidth',2);
plot(t,u1C,'b-.','LineWidth',2);
grid minor

subplot(2,1,2)
plot(t,u2A,'r:','LineWidth',2);hold on
plot(t,u2B,'g--','LineWidth',2);
plot(t,u2C,'b-.','LineWidth',2);
grid minor
legend('Nominal response without saturation','Nominal response with saturation','Response with AW compensator');
set(gcf,'Position',[23,44,1321,640]);
fCutFig(2,1)

