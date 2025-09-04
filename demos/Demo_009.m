% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        27-03-2020
% 
% -------------------------------------------------------------------------
% Demo_009:    Robust controller synthesis with LTI parametric
%              uncertainties.
%
%              This demo performs an S/KS two DoF control design for an
%              uncertain second-order system with parametric uncertainties
%              by means of the function fRobsyn.
%
% -------------------------------------------------------------------------
close all;clear all;clc;

disp('------------------------------------------------------------------');
disp('Demo_009: Robust controller synthesis with fRobsyn');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file demonstrates a Robust controller synthesis by means');
disp('of the function fRobsyn. The problem under consideration is an');
disp('S/KS 2DoF control design for an uncertain second-order system with');
disp('parametric uncertainties.');
disp(' ');
input('Press enter to continue...');
disp(' ');
disp('------------------------------------------------------------------');

% Define the uncertain plant model
k0                = 1;                                                     % Nominal time constant
d                 = 0.1;                                                   % Nominal damping
r                 = 0.4;                                                   % Variability (of time constant)
delta             = ureal('delta',0,'Range',[-1,1],'AutoSimplify','full'); % Normalized parametric uncertainty
samples           = linspace(-1,1,20);
k                 = k0*(1+r*delta);                                        % Uncertain time constant
G                 = tf(1,[k^2,2*d*k,1]);                                   % Uncertain plant model
[H,Delta]         = lftdata(G);                                            % Uncertain plant model in LFR form, G = lft(Delta,H)
H                 = ssbal(H);                                              % Balance realization matrices

% Define the open-loop system interconnection
systemnames       = 'H';
inputvar          = '[p{2};r;u]';
outputvar         = '[H(1:2);r-H(3);u;r;H(3)]';
input_to_H        = '[p(1:2);u]';
cleanupsysic      = 'yes';
olic              = sysic;

% Define weights
Wperf             = ss(-8.660258367974016e-04,1,0.865592823879067,0.5);
Wcmd              = ss(-100,32,-30.9375,10);
Wout              = blkdiag(1,1,Wperf,Wcmd,1,1);
wolic             = Wout*olic;

% Perform a nominal controller synthesis
[Knom,~,ga_nom]   = hinfsyn(wolic(3:end,3:end),2,1,'method','lmi');
disp(ga_nom)

% Define iqc uncertainty block
de                = iqcdelta('de','InputChannel',1:2,'OutputChannel',1:2,'Bounds',[-1,1]);
ude               = iqcassign(de,'ultis','Length',2);

% perform a robust controller synthesis
options.maxiter   = 5;
options.subopt    = 1.01;
options.constants = 1e-6*ones(1,3);
options.Pi11pos   = 1e-5;
options.FeasbRad  = 2e5;

[K,ga]            = fRobsyn(wolic,ude,[2,2,2],[2,1,1],options);
Krob              = K{end};

% Build closed-loop uncertain system
clic_nom          = lft(lft(Delta,olic),Knom);
clic_rob          = lft(lft(Delta,olic),Krob);
clic_nom_samp     = usubs(clic_nom,'delta',samples);
clic_rob_samp     = usubs(clic_rob,'delta',samples);

wclic_nom         = lft(lft(Delta,wolic),Knom);
wclic_rob         = lft(lft(Delta,wolic),Krob);
wclic_nom_samp    = usubs(wclic_nom,'delta',samples);
wclic_rob_samp    = usubs(wclic_rob,'delta',samples);

iter              = 1:1:length(ga);
set(0,'DefaultAxesFontSize',14);

figure(1)
bode(usubs(G,'delta',samples),'b-.',G.NominalValue,'r');
grid on
legend('Perturbed Model','Nominal Model');
title('Plant Model');
set(gcf,'Position',[23,44,1321,640]);

% Plot design weights
figure(2)
bodemag(1/Wperf,'b',1/Wcmd,'g--');
grid on
legend('1/We','1/Wu');
title('Design Weights');
set(gcf,'Position',[23,44,1321,640]);

% iteration step versus worst-case H-infinity norm gamma
figure(3)
plot(iter,ga,'or--');hold on
plot([iter(1),iter(end)],ga_nom*[1,1],'b-');hold on
xticks(iter);
axis([1,4,0.5,5.5]);
grid on
legend('K_{rob}','K_{nom}');
xlabel('Synthesis/Analysis iteration step');
ylabel('Induced L2-gain \gamma');
title('Iteration Progess');
set(gcf,'Position',[23,44,1321,640]);
fCutFig(1,1);

% Sigma plots of weighted closed-loop system
figure(4)
sigma(wclic_nom_samp,'b:',wclic_rob_samp,'r:',{1e-2,1e2});
grid on; hold all
sigma(ss(ga_nom),'b-',ss(ga(end)),'r-');
legend('Knom','Krob');
title('Sigma plots of 20 samples of the weighted closed-loop system.');
set(gcf,'Position',[23,44,1321,640]);

% Closed-loop sensitivity functions with Knom
figure(5)
subplot(2,2,1)
sigma(clic_nom_samp(1,1,:,:), 'b:', 1/Wperf, 'r--', {1e-2,1e2});
legend('Sensitivity','We^{-1}','Location','SE');
grid on
title('Sensitivity functions for 20 samples of closed-loop system with Knom.');

% Closed-loop sensitivity functions with Krob
subplot(2,2,3)
sigma(clic_rob_samp(1,1,:,:), 'b:', 1/Wperf, 'r--', {1e-2,1e2});
legend('Sensitivity','We^{-1}','Location','SE');
grid on
title('Sensitivity functions for 20 samples of closed-loop system with Krob.');

% Closed-loop KS functions with Knom
subplot(2,2,2)
sigma(clic_nom_samp(2,1,:,:), 'b:', 1/Wcmd, 'r--', {1e-2,1e3});
legend('KS','We^{-1}','Location','SE');
grid on
title('KS functions for 20 samples of closed-loop system with Knom.');

% Closed-loop sensitivity functions with Krob
subplot(2,2,4)
sigma(clic_rob_samp(2,1,:,:), 'b:', 1/Wcmd, 'r--', {1e-2,1e3});
legend('KS','We^{-1}','Location','SE');
grid on
title('KS functions for 20 samples of closed-loop system with Krob.');
set(gcf,'Position',[23,44,1321,640]);

% Step response Knom
figure(6)
subplot(1,2,1)
step(clic_nom.NominalValue(1,1),'r-',20);hold all
step(clic_nom_samp(1,1,:,:),'b:',20);
step(clic_nom.NominalValue(1,1),'r',20);
title('Reference Tracking');

% Step response Kiqc
subplot(1,2,2)
step(clic_rob.NominalValue(1,1),'r-',20);hold all
step(clic_rob_samp(1,1,:,:),'b:',20);
step(clic_rob.NominalValue(1,1),'r',20);
legend('Nominal','Peturbed');
title('Reference Tracking');
set(gcf,'Position',[23,44,1321,640]);

% Remove fast poles/zeros
Krob_r            = fRFZP(Krob,3e3);
clic_rob_r        = lft(lft(Delta,olic),Krob_r);
clic_rob_r_samp   = usubs(clic_rob_r,'delta',samples);
wclic_rob_r       = lft(lft(Delta,wolic),Krob_r);
[wcg_r,wcu_r]     = wcgain(wclic_rob_r);

% Step response Knom
figure(7)
subplot(1,2,1)
step(clic_nom.NominalValue(1,1),'r-',20);hold all
step(clic_nom_samp(1,1,:,:),'b:',20);
step(clic_nom.NominalValue(1,1),'r',20);
title('Reference Tracking');

% Step response Kiqc
subplot(1,2,2)
step(clic_rob_r.NominalValue(1,1),'r-',20);hold all
step(clic_rob_r_samp(1,1,:,:),'b:',20);
step(clic_rob_r.NominalValue(1,1),'r',20);
legend('Nominal','Peturbed');
title('Reference Tracking');
set(gcf,'Position',[23,44,1321,640]);