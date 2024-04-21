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
% Date:        20-03-2020
% 
% -------------------------------------------------------------------------
% Demo_007:    H-infinity controller synthesis with fHinfsyn
%
%              This demo demonstrates and validates the function "fHinfsyn"
%              with the function "hinfsyn" by performing a controller
%              synthesis for the active suspension demo available in
%              MATLAB.
% -------------------------------------------------------------------------
close all;clear all;clc;set(0,'DefaultAxesFontSize',14);

clc
disp('------------------------------------------------------------------');
disp('Demo_007: H-infinity controller synthesis with fHinfsyn');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file demonstrates an H-infinity synthesis by means of');
disp('the function fHinfsyn. The results are compared with with the');
disp('H-infinity and mu-synthesis algorithms that are available in the');
disp('Robust Control Toolbox.');
disp(' ');
input('Press enter to continue...');
disp(' ');
disp('------------------------------------------------------------------');

% Physical parameters
mb                         = 300;    % kg
mw                         = 60;     % kg
bs                         = 1000;   % N/m/s
ks                         = 16000 ; % N/m
kt                         = 190000; % N/m

% State space matrices
A                          = [0,1,0,0;[-ks,-bs,ks,bs]/mb;0,0,0,1;[ks,bs,-ks-kt,-bs]/mw];
B                          = [0,0;0,1e3/mb;0,0;[kt,-1e3]/mw];
C                          = [1,0,0,0;1,0,-1,0;A(2,:)];
D                          = [0,0;0,0;B(2,:)];
qcar                       = ss(A,B,C,D);
qcar.StateName             = {'body travel (m)';'body vel (m/s)';'wheel travel (m)';'wheel vel (m/s)'};
qcar.InputName             = {'r';'fs'};
qcar.OutputName            = {'xb';'sd';'ab'};

figure(1)
bodemag(qcar({'ab','sd'},'r'),'b',qcar({'ab','sd'},'fs'),'r--',{1 100});
legend('Road disturbance (r)','Actuator force (fs)','location','SouthWest')
title({'Gain from road dist (r) and actuator force (fs) to body accel (ab) and suspension travel (sd)'});
set(gcf,'Position',[23,44,1321,640]);
   
ActNom                     = tf(1,[1/60,1]);
Wunc                       = makeweight(0.40,15,3);
unc                        = ultidyn('unc',[1 1],'SampleStateDim',5);
Act                        = ActNom*(1 + Wunc*unc);
Act.InputName              = 'u';
Act.OutputName             = 'fs';
rng('default')

figure(2)
bode(Act,'b:',Act.NominalValue,'r-',logspace(-1,3,120));
grid on
set(gcf,'Position',[23,44,1321,640]);

Wroad                      = ss(0.07);
Wroad.u                    = 'd1';
Wroad.y                    = 'r';

Wact                       = 0.8*tf([1 50],[1 500]);
Wact.u                     = 'u';
Wact.y                     = 'e1';

Wd2                        = ss(0.01);
Wd2.u                      = 'd2';
Wd2.y                      = 'Wd2';

Wd3                        = ss(0.5);
Wd3.u                      = 'd3';
Wd3.y                      = 'Wd3';

HandlingTarget             = 0.04 * tf([1/8 1],[1/80 1]);
ComfortTarget              = 0.4 * tf([1/0.45 1],[1/150 1]);
Targets                    = [HandlingTarget;ComfortTarget];

figure(3)
bodemag(qcar({'sd','ab'},'r')*Wroad,'b',Targets,'r--',{1,1000});
grid on
title('Response to road disturbance');
legend('Open-loop','Closed-loop target');
set(gcf,'Position',[23,44,1321,640]);

% Three design points
beta                       = reshape([0.01 0.5 0.99],[1 1 3]);

Wsd                        = beta/HandlingTarget;
Wsd.u                      = 'sd';
Wsd.y                      = 'e3';

Wab                        = (1-beta)/ComfortTarget;
Wab.u                      = 'ab';
Wab.y                      = 'e2';

sdmeas                     = sumblk('y1 = sd+Wd2');
abmeas                     = sumblk('y2 = ab+Wd3');

ICinputs                   = {'d1';'d2';'d3';'u'};
ICoutputs                  = {'e1';'e2';'e3';'y1';'y2'};

qcaric                     = connect(qcar(2:3,:),Act,Wroad,Wact,Wab,Wsd,Wd2,Wd3,sdmeas,abmeas,ICinputs,ICoutputs);
             
ncont                      = 1; % one control signal, u
nmeas                      = 2; % two measurement signals, sd and ab
K1                         = ss(zeros(ncont,nmeas,3));
K2                         = ss(zeros(ncont,nmeas,3));
gamma1                     = zeros(3,1);
gamma2                     = zeros(3,1);

for i = 1:3
   [K1(:,:,i),~,gamma1(i)] = hinfsyn(minreal(qcaric(:,:,i).NominalValue),nmeas,ncont);
   options.subopt          = 1.02;
   options.condnr          = 1.01;
   options.FeasbRad        = 1e6;
   [K2(:,:,i),gamma2(i)]   = fHinfsyn(minreal(qcaric(:,:,i).NominalValue),nmeas,ncont,options);
end

disp('Hinf-norm weighted closed-loop system obtained by "hinfsyn":');
disp(gamma1);
disp('Hinf-norm weighted closed-loop system obtained by "fHinfsyn":');
disp(gamma2);

figure(4)
subplot(211)
bodemag(K1(:,:,1),'b-',K1(:,:,2),'r--',K1(:,:,3),'g-.');
title('Bode magnitude plot of K (designed by hinfsyn)');
grid on
subplot(212)
bodemag(K2(:,:,1),'b-',K2(:,:,2),'r--',K2(:,:,3),'g-.');
title('Bode magnitude plot of K (designed by fHinfsyn)');
grid on
legend('Comfort','Balanced','Handling','location','SouthEast')
set(gcf,'Position',[23,44,1321,640]);

% Closed-loop models
K1.u                       = {'sd','ab'};
K1.y                       = 'u';
CL1                        = connect(qcar,Act.Nominal,K1,'r',{'xb';'sd';'ab'});

K2.u                       = {'sd','ab'};
K2.y                       = 'u';
CL2                        = connect(qcar,Act.Nominal,K2,'r',{'xb';'sd';'ab'});

% Plot results for K1
figure(5)
subplot(3,2,1)
bodemag(qcar('xb','r'),'b',CL1(1,1,1),'r-.',CL1(1,1,2),'m-.',CL1(1,1,3),'k-.',{1,140});
title('Body travel due to road (for hinfsyn)');
grid on

subplot(3,2,2)
bodemag(qcar('xb','r'),'b',CL2(1,1,1),'r-.',CL2(1,1,2),'m-.',CL2(1,1,3),'k-.',{1,140});
title('Body travel due to road (for fHinfsyn)');
grid on

subplot(3,2,3)
bodemag(qcar('sd','r'),'b',CL1(2,1,1),'r-.',CL1(2,1,2),'m-.',CL1(2,1,3),'k-.',{1,140});
title('Suspension deflection due to road (for hinfsyn)');
grid on

subplot(3,2,4)
bodemag(qcar('sd','r'),'b',CL2(2,1,1),'r-.',CL2(2,1,2),'m-.',CL2(2,1,3),'k-.',{1,140});
title('Suspension deflection due to road (for fHinfsyn)');
grid on

subplot(3,2,5)
bodemag(qcar('ab','r'),'b',CL1(3,1,1),'r-.',CL1(3,1,2),'m-.',CL1(3,1,3),'k-.',{1,140});
title('Body acceleration due to road (for hinfsyn)');
grid on

subplot(3,2,6)
bodemag(qcar('ab','r'),'b',CL2(3,1,1),'r-.',CL2(3,1,2),'m-.',CL2(3,1,3),'k-.',{1,140});
title('Body acceleration due to road (for fHinfsyn)');
grid on
legend('Open-loop','Comfort','Balanced','Handling','location','SouthEast')
set(gcf,'Position',[23,44,1321,800]);

% Road disturbance
t                          = 0:0.0025:1;
roaddist                   = zeros(size(t));
roaddist(1:101)            = 0.025*(1-cos(8*pi*t(1:101)));

% Closed-loop model
SIMK1                      = connect(qcar,Act.Nominal,K1,'r',{'xb';'sd';'ab';'fs'});
SIMK2                      = connect(qcar,Act.Nominal,K2,'r',{'xb';'sd';'ab';'fs'});

% Simulate
p11                        = lsim(qcar(:,1),roaddist,t);
y11                        = lsim(SIMK1(1:4,1,1),roaddist,t);
y21                        = lsim(SIMK1(1:4,1,2),roaddist,t);
y31                        = lsim(SIMK1(1:4,1,3),roaddist,t);

p12                        = lsim(qcar(:,1),roaddist,t);
y12                        = lsim(SIMK2(1:4,1,1),roaddist,t);
y22                        = lsim(SIMK2(1:4,1,2),roaddist,t);
y32                        = lsim(SIMK2(1:4,1,3),roaddist,t);

% Plot results for K1
figure(6)
subplot(221)
plot(t,p11(:,1),'b',t,y11(:,1),'r.',t,y21(:,1),'m.',t,y31(:,1),'k.',t,roaddist,'g');
title('Body travel'), ylabel('x_b (m)');

subplot(222)
plot(t,p11(:,2),'b',t,y11(:,2),'r.',t,y21(:,2),'m.',t,y31(:,2),'k.',t,roaddist,'g');
title('Suspension deflection'), xlabel('Time (s)'), ylabel('s_d (m)');

subplot(223)
plot(t,p11(:,3),'b',t,y11(:,3),'r.',t,y21(:,3),'m.',t,y31(:,3),'k.',t,roaddist,'g');
title('Body acceleration'), ylabel('a_b (m/s^2)');

subplot(224)
plot(t,zeros(size(t)),'b',t,y11(:,4),'r.',t,y21(:,4),'m.',t,y31(:,4),'k.',t,roaddist,'g');
title('Control force'), xlabel('Time (s)'), ylabel('f_s (kN)');
legend('Open-loop','Comfort','Balanced','Handling','Road Disturbance','location','SouthEast')
set(gcf,'Position',[23,44,1321,640]);

% Plot results for K2
figure(7)
subplot(221)
plot(t,p12(:,1),'b',t,y12(:,1),'r.',t,y22(:,1),'m.',t,y32(:,1),'k.',t,roaddist,'g');
title('Body travel'), ylabel('x_b (m)');

subplot(222)
plot(t,p12(:,2),'b',t,y12(:,2),'r.',t,y22(:,2),'m.',t,y32(:,2),'k.',t,roaddist,'g');
title('Suspension deflection'), xlabel('Time (s)'), ylabel('s_d (m)');


subplot(223)
plot(t,p12(:,3),'b',t,y12(:,3),'r.',t,y22(:,3),'m.',t,y32(:,3),'k.',t,roaddist,'g');
title('Body acceleration'), ylabel('a_b (m/s^2)');

subplot(224)
plot(t,zeros(size(t)),'b',t,y12(:,4),'r.',t,y22(:,4),'m.',t,y32(:,4),'k.',t,roaddist,'g');
title('Control force'), xlabel('Time (s)'), ylabel('f_s (kN)');
legend('Open-loop','Comfort','Balanced','Handling','Road Disturbance','location','SouthEast');
set(gcf,'Position',[23,44,1321,640]);

[Krob1,rpMU]               = musyn(qcaric(:,:,2),nmeas,ncont);
Krob1.u                    = {'sd','ab'};
Krob1.y                    = 'u';

[P,Delta]                  = lftdata(qcaric(:,:,2));
options.subopt             = 1.02;
options.condnr             = 1.01;
options.FeasbRad           = 1e6;
options.alpha              = 1;
[Krob2,ga2]                = fHinfsyn(P,[1,3,2],[1,3,1],options);
Krob2.u                    = {'sd','ab'};
Krob2.y                    = 'u';
[wcg_nom,wcu_nom]          = wcgain(lft(qcaric(:,:,2),K1));
[wcg_rob1,wcu_rob1]        = wcgain(lft(qcaric(:,:,2),Krob1));
[wcg_rob2,wcu_rob2]        = wcgain(lft(qcaric(:,:,2),Krob2));

% Closed-loop model (nominal)
SIMKrob1                   = connect(qcar,Act.Nominal,Krob1,'r',{'xb';'sd';'ab';'fs'});
SIMKrob2                   = connect(qcar,Act.Nominal,Krob2,'r',{'xb';'sd';'ab';'fs'});

% Simulate
p1                         = lsim(qcar(:,1),roaddist,t);
y1                         = lsim(SIMKrob1(1:4,1),roaddist,t);
y2                         = lsim(SIMKrob2(1:4,1),roaddist,t);

% Plot results
figure(8)
subplot(221)
plot(t,p1(:,1),'b',t,y1(:,1),'r',t,roaddist,'g');
title('Body travel'), ylabel('x_b (m)')

subplot(222)
plot(t,p1(:,3),'b',t,y1(:,3),'r')
title('Body acceleration'), ylabel('a_b (m/s^2)')

subplot(223)
plot(t,p1(:,2),'b',t,y1(:,2),'r')
title('Suspension deflection'), xlabel('Time (s)'), ylabel('s_d (m)')

subplot(224)
plot(t,zeros(size(t)),'b',t,y1(:,4),'r')
title('Control force'), xlabel('Time (s)'), ylabel('f_s (kN)')
legend('Open-loop','Robust design','location','SouthEast')
set(gcf,'Position',[23,44,1321,640]);

% Plot results
figure(9)
subplot(221)
plot(t,p1(:,1),'b',t,y2(:,1),'r',t,roaddist,'g');
title('Body travel'), ylabel('x_b (m)')

subplot(222)
plot(t,p1(:,3),'b',t,y2(:,3),'r')
title('Body acceleration'), ylabel('a_b (m/s^2)')

subplot(223)
plot(t,p1(:,2),'b',t,y2(:,2),'r')
title('Suspension deflection'), xlabel('Time (s)'), ylabel('s_d (m)')

subplot(224)
plot(t,zeros(size(t)),'b',t,y2(:,4),'r')
title('Control force'), xlabel('Time (s)'), ylabel('f_s (kN)')
legend('Open-loop','Robust design','location','SouthEast')
set(gcf,'Position',[23,44,1321,640]);

rng('default');
nsamp                      = 100;
CLU1                       = connect(qcar,Act,K1(:,:,2),'r',{'xb','sd','ab'});
CLUrob1                    = connect(qcar,Act,Krob1,'r',{'xb','sd','ab'});
CLUrob2                    = connect(qcar,Act,Krob2,'r',{'xb','sd','ab'});

% Uncertain closed-loop model with balanced H-infinity controller
figure(10)
lsim(usample(CLU1,nsamp),'b',CLU1.Nominal,'r',roaddist,t)
title('Nominal "balanced" design: K1')
legend('Perturbed','Nominal','location','SouthEast')
set(gcf,'Position',[23,44,1321,640]);

% Uncertain closed-loop model with balanced robust controller
figure(11)
lsim(usample(CLUrob1,nsamp),'b',CLUrob1.Nominal,'r',roaddist,t)
title('Robust "balanced" design: Krob1');
legend('Perturbed','Nominal','location','SouthEast');
set(gcf,'Position',[23,44,1321,640]);

% Uncertain closed-loop model with balanced robust controller
figure(12)
lsim(usample(CLUrob2,nsamp),'b',CLUrob2.Nominal,'r',roaddist,t);
title('Robust "balanced" design: Krob2');
legend('Perturbed','Nominal','location','SouthEast');
set(gcf,'Position',[23,44,1321,640]);