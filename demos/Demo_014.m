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
% Date:        25-08-2023
% 
% -------------------------------------------------------------------------
% Demo_014:    Robust estimator design
% 
%              This demo demonstrates how to design a robust estimator
%              using fRobest
% -------------------------------------------------------------------------
clc;close all;

%% Define plant
A                 = [0,0,1,0;0,0,0,1;-2,1,-2,2;2,-2,4,-4];
Bp                = [zeros(2);-1.75,-0.5;3.5,1];
Bw                = [zeros(2);1,0;0,0];
Cq                = [0,0,1,-1;1,-1,0,0];
Dqp               = zeros(2);
Dqw               = zeros(2);
Cz                = [0,1,0,0];
Dzp               = zeros(1,2);
Dzw               = zeros(1,2);
Cy                = [1,0,0,0];
Dyp               = zeros(1,2);
Dyw               = [0,1];
H                 = ss(A,[Bp,Bw],[Cq;Cz;Cy],[Dqp,Dqw;Dzp,Dzw;Dyp,Dyw]);

%% Nominal H2 - Estimator design
systemnames       = 'H';
inputvar          = '[p{2};w{2};u]';
outputvar         = '[H(1:2);H(3);u;H(3)+u;H(4)]';
input_to_H        = '[p;w]';
cleanupsysic      = 'yes';
olic              = sysic;

[E,~,gamnom]      = h2syn(olic(5:end,3:end),1,1);

disp('Computed H2-norm:');
disp(gamnom)

Enom              = balreal(E);

%% Robustness H2 analysis for Enom
CL                = balreal(lft(olic([1:2,5:6],:),Enom));
al                = linspace(0,1,21);al(1) = 1e-4;

for i = 1:length(al)
    de1           = iqcdelta('de1','InputChannel',1,'OutputChannel',1,'Bounds',al(i)*[-1,1]);
    de1           = iqcassign(de1,'ultis','Length',3);

    de2           = iqcdelta('de2','InputChannel',2,'OutputChannel',2,'Bounds',al(i)*[-1,1]);
    de2           = iqcassign(de2,'ultis','Length',3);

    pe            = iqcdelta('pe','ChannelClass','P','InputChannel',3:4,'OutputChannel',3,'PerfMetric','H2');

    Delta         = {de1,de2,pe};
    prob          = iqcanalysis(CL,Delta);

    gamnom(i)     = prob.gamma;
end
figure(3)
plot(al,gamnom,'r');hold on

%% Robust H2 estimator design
de1               = iqcdelta('de1','InputChannel',1,'OutputChannel',1,'Bounds',[-1,1]);
de1               = iqcassign(de1,'ultis','Length',3);

de2               = iqcdelta('de2','InputChannel',2,'OutputChannel',2,'Bounds',[-1,1]);
de2               = iqcassign(de2,'ultis','Length',3);

options.perf      = 'H2';
options.StrProp   = 'yes';
options.subopt    = 1.03;
options.constants = 1e-8*ones(1,4);
options.Pi11pos   = 1e-8;
options.FeasbRad  = 1e9;

[Erob,gamrob]     = fRobest(H,{de1,de2},[2,1,1],[2,2],options);

disp('Computed worst case H2-norm:');
disp(gamrob)

om                = logspace(-2,2,500);
[ma_nom,~]        = bode(Enom,om);
[ma_rob,~]        = bode(Erob,om);


figure(4)
semilogx(om,20*log10(squeeze(ma_nom)));hold on
semilogx(om,20*log10(squeeze(ma_rob)));hold on
title('Bode magnitude plots of E_{nom} and E_{rob}');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on
legend('Enom','Erob');
set(gcf,'Position',[300,400,700,350]);
fCutFig(1,1);

%% Robustness H2 analysis for Erob
CL                = balreal(lft(olic([1:2,5:6],:),Erob));
al                = linspace(0,1,21);al(1) = 1e-4;

for i = 1:length(al)
    de1           = iqcdelta('de1','InputChannel',1,'OutputChannel',1,'Bounds',al(i)*[-1,1]);
    de1           = iqcassign(de1,'ultis','Length',3);

    de2           = iqcdelta('de2','InputChannel',2,'OutputChannel',2,'Bounds',al(i)*[-1,1]);
    de2           = iqcassign(de2,'ultis','Length',3);
    
    pe            = iqcdelta('pe','ChannelClass','P','InputChannel',3:4,'OutputChannel',3,'PerfMetric','H2');

    Delta         = {de1,de2,pe};
    prob          = iqcanalysis(CL,Delta);

    gamrob(i)     = prob.gamma;
end

figure(3)
plot(al,gamrob,'k--');
xlabel('\alpha');
ylabel('\gamma');
title('Computed worst case H2-norm \gamma for increasing values of |\delta_i|<\alpha, i=1,2');
legend('E_{nom}','E_{rob}');
set(gcf,'Position',[300,400,700,350]);
fCutFig(1,1);

%% Time domain simulations
try;close(5);catch;end
try;close(6);catch;end

rng(132);
tend              = 1e4;
t                 = linspace(0,tend,1e5);
u1                = ureal('u1',0);
u2                = ureal('u2',0);
ax1               = 0.25;
ax2               = 0.75;

CLnom             = lft(lft(blkdiag(u1,u2),olic),Enom);
CLrob             = lft(lft(blkdiag(u1,u2),olic),Erob);

unc               = linspace(-1,1,5);

for i = 1:5
    CLnom_sub     = usubs(CLnom,'u1',unc(i),'u2',unc(i));
    CLrob_sub     = usubs(CLrob,'u1',unc(i),'u2',unc(i));
    
    w1            = 0.01*randn(1,length(t));

    s             = tf('s');
    W             = 1/(s/0.0001+1);
    w1f           = 1e4*lsim(W,w1,t);

    w2            = 0.00001*randn(1,length(t));
    enom          = lsim(CLnom_sub,[w1f';w2],t);
    erob          = lsim(CLrob_sub,[w1f';w2],t);
    
    RMSquare(1,i) = rms(enom(:,3));
    RMSquare(2,i) = rms(erob(:,3));

    figure(5)
    subplot(2,5,1:5);
    plot(t,enom(:,3));hold on
    
    subplot(2,5,5+i);
    plot(t,enom(:,1),'k--');hold on
    plot(t,-enom(:,2),'r-');

    if i == 5
        subplot(2,5,1:5);
        title('Estimation errors E_{nom} for different values of \delta_1=\delta_2')
        axis([0,tend,-ax1,ax1])
        xlabel('Time (s)')
        ylabel('e')

        subplot(2,5,6);
        ylabel('z and zhat');
        axis([0,tend,-ax2,ax2])

        subplot(2,5,7);
        axis([0,tend,-ax2,ax2])

        subplot(2,5,8);
        xlabel('Time (s)')
        axis([0,tend,-ax2,ax2])

        subplot(2,5,9);
        axis([0,tend,-ax2,ax2])

        subplot(2,5,10);
        axis([0,tend,-ax2,ax2])
    end

    figure(6)
    subplot(2,5,1:5);
    plot(t,erob(:,3));hold on
    
    subplot(2,5,5+i);
    plot(t,erob(:,1),'k--');hold on
    plot(t,-erob(:,2),'r-');

    if i == 5
        subplot(2,5,1:5);
        axis([0,tend,-ax1,ax1])
        title('Estimation errors E_{rob} for different values of \delta_1=\delta_2')
        xlabel('Time (s)')
        ylabel('e')

        subplot(2,5,6);
        ylabel('z and zhat');
        axis([0,tend,-ax2,ax2])

        subplot(2,5,7);
        axis([0,tend,-ax2,ax2])

        subplot(2,5,8);
        xlabel('Time (s)')
        axis([0,tend,-ax2,ax2])

        subplot(2,5,9);
        axis([0,tend,-ax2,ax2])

        subplot(2,5,10);
        axis([0,tend,-ax2,ax2])
    end
end

disp('RMS values for the performed simulations');
disp('Nominal estimator design Enom')
disp(RMSquare(1,:))

disp('Robust estimator design Erob')
disp(RMSquare(2,:))

disp('Mean RMS values for the performed simulations');
disp('Nominal estimator design Enom')
disp(mean(RMSquare(1,:)))

disp('Robust estimator design Erob')
disp(mean(RMSquare(2,:)))

figure(5)
set(gcf,'Position',[300,400,700,350]);

figure(6)
set(gcf,'Position',[300,400,700,350]);
