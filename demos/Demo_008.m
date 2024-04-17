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
% Date:        24-03-2020
% 
% -------------------------------------------------------------------------
% Demo_008:    LPV controller synthesis with fLPVsyn
%
%              This demo performs an LPV controller synthesis for a missile
%              control design problem.
%
% -------------------------------------------------------------------------
close all;clear all;clc;set(0,'DefaultAxesFontSize',14);

disp('------------------------------------------------------------------');
disp('Demo_008: LPV controller synthesis with fLPVsyn');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file demonstrates an LPV controller synthesis by means');
disp('of the function fLPVsyn. The problem under consideration is an');
disp('parametric flight control problem with the Mach number and flight');
disp('path angle as scheduling paramters.');
disp(' ');
input('Press enter to continue...');
disp(' ');
disp('------------------------------------------------------------------');

% Load missile open-loop plant
load missile

% Define uncertainties
del_a                      = ureal('del_a',0,'Range',[-1,1]);
del_M                      = ureal('del_M',0,'Range',[-1,1]);
Del_a                      = pi/18*(1 + del_a);
Del_M                      = 3 + del_M;
Del_T                      = blkdiag(Del_a*eye(2),Del_M*eye(5));

% Extract open-loop interconnection (olic) for normalized uncertainties
P                          = lft(Del_T,ssbal(G));
[olic,De]                  = lftdata(P);

% Append weighting functions
wp                         = ss(tf([0.5,3.5],[1,3.5e-6]));
wu                         = ss(20*tf([1/0.9,1],[1e-6,1]));
Wo                         = ssbal(blkdiag(eye(7),wp,wu,eye(2)));
wolic                      = Wo*ssbal(olic);

% Define scheduling block
H{1}                       = blkdiag(eye(5),zeros(2));
H{2}                       = blkdiag(zeros(5),eye(2));
La                         = polydec(pvec('box',[-1,1;-1,1]))';
iqcdeltaOpt.InputChannel   = 1:7;
iqcdeltaOpt.OutputChannel  = 1:7;
iqcdeltaOpt.Polytope       = La;
iqcdeltaOpt.UncertaintyMap = H;
iqcdeltaOpt.TimeInvTimeVar = 'TV';
iqcdeltaOpt.Structure      = 'FB';
Delta                      = iqcdelta('Delta',iqcdeltaOpt);

% Define synthesis options
options.RelaxType          = 'CH';
options.FeasbRad           = 1e8;
options.constants          = [1e-8,1e-7,1e-7,1e-9,1e-8,1e-8];
options.subopt             = 1.1;
options.bounds             = [1e5,1,1e5,1];

% Perform synthesis
[K,ga]                     = fLPVsyn(wolic,Delta,[7,2,2],[7,1,1],options);

% Construct closed-loop plant
De_S                       = fOblkdiag(De');
Ksch                       = lft(K,De_S);
Pcl                        = lft(P,Ksch);

% simulate results
sim('missile');

figure(1)
subplot(2,1,1)
plot(Nz.Time,Nz.Data(:,1),'k:','LineWidth',1);hold on
plot(Nz.Time,Nz.Data(:,2),'r-','LineWidth',2);hold on
set(gcf,'Position',[23,44,1321,640]);
xlabel('Time (s)');
ylabel('Normal acc/g (n_z)');

subplot(2,1,2)
plot(Mach.Time,Mach.Data(:,1),'b-','LineWidth',1);
ylabel('Mach nr.');
set(gcf,'Position',[23,44,1321,640]);
fCutFig(2,1);
