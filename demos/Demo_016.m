% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        30-08-2025
% 
% -------------------------------------------------------------------------
% Demo_016:    Data-based parameter refinement
%
%              This demo file performs a data-based parameter refinement
%              for a controlled uncertain system with two uncertainties.
%              The example proceeds in 3 steps:
%                1.) Perform a nominal controller synthesis. This
%                    controllwe works well for the nominal scenario, while
%                    the performance degrades for various combinations of
%                    the uncertainties.
%                2.) Collect input-output data from an experiment and run
%                    the function fRefine with the aim to reduce the
%                    uncertainty bounds as much as possible
%                3.) Based on the reduced parameter bounds, compute the new
%                    nominal (mean value of the uncertainties) system and
%                    re-run the synthesis. The new controller will be
%                    insensitive to variations of the uncertainties due to
%                    the (much) reduced parameter bounds.
%
%              As input, it is possible to specify the following:
%
%                1.) Both parameters are to be refined'
%                2.) Only the parameter k is to be refined while d is
%                    assumed to be nominal.
%                3.) Only the parameter k is to be refined while d is
%                    assumed to be uncertain.
%
% -------------------------------------------------------------------------
close all;clc;set(0,'DefaultAxesFontSize',14);
disp('------------------------------------------------------------------');
disp('Demo_016: Parametric sensitivity analysis                         ');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file performs a data-based parameter refinement for a   ');
disp('controlled uncertain system with two uncertainties. See the       ');
disp('description at the beginning of this file for further details.    ');
disp(' ');
disp('------------------------------------------------------------------');
disp('Please specify the configuration: ');
disp('  1: Both parameters are to be refined.');
disp('  2: Only the parameter k is to be refined while d is assume to be nominal.');
disp('  3: Only the parameter k is to be refined while d is assumed to be uncertain.');
disp(' ');
me = input('Choose option 1, 2 or 3: ');

% Define the uncertain plant model
d0                     = 0.1;                                   % Nominal damping
dd                     = 0.4;
d                      = d0*(1+dd*ureal('d',0,'Range',[-1,1],...
                        'AutoSimplify','full'));                % Normalized parametric uncertainty
k0                     = 1;                                     % Nominal time constant
kd                     = 0.4;                                   % Variability (of time constant)
k                      = k0*(1+kd*ureal('k',0,'Range',[-1,1],...
                        'AutoSimplify','full'));                % Normalized parametric uncertainty
samples                = linspace(-1,1,20);
G                      = tf(1,[k^2,2*d*k,1]);                   % Uncertain plant model
[H,De]                 = lftdata(G);                            % Uncertain plant model in LFR form, G = lft(Delta,H)
H                      = ssbal(H);                              % Balance realization matrices

% Define the open-loop system interconnection
systemnames            = 'H';
inputvar               = '[p{4};n;r;u]';
outputvar              = '[H(1:4);r-H(5);u;H(5)+n;r;H(5)+n]';
input_to_H             = '[p(1:4);u]';
cleanupsysic           = 'yes';
P                      = sysic;

% Define weights for controller synthesis
We                     = ss(-8.660258367974016e-04,1,0.865592823879067,0.5);
Wu                     = ss(-100,32,-30.9375,10);
% Wn                     = ss(-1,0.03125,-0.032,0.001);
Wn                     = ss(-1,0.125,-0.08,0.01);
Wout                   = blkdiag(1,1,1,1,We,Wu,1,1,1);
Win                    = blkdiag(1,1,1,1,Wn,1,1);
Pw                     = Wout*P*Win;

% Perform a nominal controller synthesis
[Knom,~,ga_nom]        = hinfsyn(Pw([5,6,8,9],5:end),2,1,'method','lmi');
disp(['gamma nominal controller = ',num2str(ga_nom)]);disp(' ');

% Postpro nominal controller by removing poles/zeros at higher frequencies
Knomred                = fRFZP(Knom,25);

% Discretize closed-loop plant
N                      = lft(lft(De,P),Knomred);
[M,De]                 = lftdata(N);
Md                     = c2d(M,0.1);
Nd                     = lft(De,Md);

k_true                 = 0.635;
if me == 2
    d_true             = 0;
else
    d_true             = -0.398;
end

Nd_true                = usubs(Nd,'k',k_true,'d',d_true);

rng(2, "twister");                           % For reproducibility
t                      = 0:0.1:20;           % Time horizon
hor                    = length(t);          % Horizon length
eps                    = 0.001;              % Process noise magnitude
n                      = eps * rand(1, hor); % Process noise
h6                     = floor(hor/6);
r                      = 5*[0.5*randn(1, h6), -4*ones(1, h6), 2*randn(1, h6), ...
                         2*ones(1, h6), -3*randn(1, h6), rand(1, h6+3)]; % Some input

% Simulation of this system
[y, ~, x]              = lsim(Nd_true,[n;r],t,zeros(size(Nd_true.a,1),1));

sigs.t                 = t;
sigs.xt                = x;
sigs.rt                = r';
sigs.yt                = y;
sigs.nt                = eps;

%%
Delreal                = [k_true,d_true];
if me == 1
    params             = {'k','d'};
elseif me == 2
    params             = {'k'};
    Nd                 = usubs(Nd,'d',0);
elseif me == 3
    params             = {'k'};
end
options.Case           = 'IOdata';
options.LiftingNo      = 3;

tic
[Nr,Bnds]              = fRefine(Nd,sigs,params,options);
disp(' ');
toc

close all
lb = length(Bnds);
clear bnds
for j = 1:size(Bnds{1},1)
    for i = 1:lb
        bnds{j}(i,1:2) = Bnds{i}(j,1:2);
    end
    figure(j)
    plot([sigs.t(1),sigs.t(end)],Delreal(j)*[1,1],'k--');hold on
    plot([sigs.t(1),sigs.t(end)],[1,1],'k-','LineWidth',2);hold on
    plot(sigs.t(1:lb),bnds{j},'r--','LineWidth',2);
    plot([sigs.t(1),sigs.t(end)],[-1,-1],'k-','LineWidth',2);
    grid on
    ylim([-1.1,1.1]);
    xlabel('Time (s)');
    title(params{j});
    legend('True value','Initial bounds','Refined bounds');
    set(gcf,'Position',[300,400,700,350]);
end

% Perform controller synthesis for refined uncertain system
if abs(Bnds{end}(1)-Bnds{end}(2)) < 1e-6
   k_red               = 0; 
else
    k_red              = ureal('k_red',mean(Bnds{end}(1,:)),'Range',Bnds{end}(1,:));
end
if me == 1
    d_red              = ureal('d_red',mean(Bnds{end}(2,:)),'Range',Bnds{end}(2,:));
elseif me == 2
    d_red              = 0;
elseif me == 3
    d_red              = d;
end
Pred                   = lft(blkdiag(d_red,k_red*eye(3)),P);
Predw                  = Wout(5:end,5:end)*Pred*Win(5:end,5:end);

% Perform a nominal controller synthesis
[Kref,~,ga_ref]        = hinfsyn(Predw.NominalValue([1,2,4,5],:),2,1,'method','lmi');
disp(' ');disp(['gamma retuned controller = ',num2str(ga_ref)]);disp(' ');

Krefred                = fRFZP(Kref,25);
Nref                   = lft(Pred,Krefred);

% Plot step responses
Ns                     = 20;
Nsamp                  = usample(N(1,2),Ns);
Nsamp(:,:,1,1)         = N.NominalValue(1,2);
Nrefsamp               = usample(Nref(1,2),Ns);
Nrefsamp(:,:,1,1)      = Nref.NominalValue(1,2);

f1                     = j+1;
f2                     = j+2;

for i = 1:Ns
    yN                 = step(Nsamp(:,:,i,1),t);
    yNref             = step(Nrefsamp(:,:,i,1),t);
    
    if i == 1
        figure(f1);
        plot(t(:),yN,'r','LineWidth',2);hold on
        figure(f2);
        plot(t(:),yNref,'r','LineWidth',2);hold on
    else
        figure(f1);
        plot(t(:),yN,'b');hold on
        figure(f2);
        plot(t(:),yNref,'b');hold on
    end
end
figure(f1)
xlabel('Time (s)');
xlim([0,10]);
title('Tracking error step responses original system');
legend('Nominal response','uncertain responses');
grid on
set(gcf,'Position',[300,400,700,350]);

figure(f2)
xlabel('Time (s)');
xlim([0,10]);
title('Tracking error step responses refined system');
legend('Nominal response','uncertain responses');
grid on
set(gcf,'Position',[300,400,700,350]);
