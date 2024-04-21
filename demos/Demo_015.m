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
% Date:        10-04-2024
% 
% -------------------------------------------------------------------------
% Demo_015:    Sensitivity analysis for uncertain LTI systems
%
%              This demo file performs a sensitivity analysis based on
%              Morris or the variance-based method for a satellite attitude
%              control problem with various uncertainties (Inertia,
%              Reaction wheel (RWL) scale factor, RWL misalignments).
% 
%              Exploiting the outcome of the sensitivity analysis the code
%              compares the following upper-bound analyses:
%
%                - Nominal performance (no uncertainty)
%                - Worst-case performance for the subset of parameters that
%                  do not affect the performance of the system (this should
%                  yield a value that is close to the nominal performance
%                  level) 
%                - Worst-case performance for the subset of parameters that
%                  do affect the performance of the system (this should
%                  yield a value that is close to the worst-case
%                  performance level)
%                - Worst-case performance for the full set of
%                  uncertainties.
%
%              As input, it is possible to specify the following mehtods:
%
%                - 1 for 'morris'
%                - 2 for 'var'
%
% -------------------------------------------------------------------------
close all;clc;
disp('------------------------------------------------------------------');
disp('Demo_015: Parametric sensitivity analysis                         ');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo file performs a sensitivity analysis based on Morris or ');
disp('the variance-based method for a satellite attitude control problem');
disp('with various uncertainties. Here it is possible to specify the    ');
disp('analysis method.                                                  ');
disp(' ');
disp('------------------------------------------------------------------');
disp('Please specify the analysis method:');
disp('  1: Morris method');
disp('  2: Variance based-approach');
disp(' ');
me = input('Choose 1 or 2: ');

if me == 1
    meth = 'morris';
elseif me == 2
    meth = 'var';
else
    meth = 'morris';
end

% Define plant model
s               = tf('s');

wd11            = 0.01;                                                                 % uncertainty weight inetria term J11
wd12            = 0.5;                                                                  % uncertainty weight inetria term J12
wd13            = 0.01;                                                                 % uncertainty weight inetria term J13
wd22            = 0.1;                                                                  % uncertainty weight inetria term J22
wd23            = 0.01;                                                                 % uncertainty weight inetria term J23
wd33            = 0.5;                                                                  % uncertainty weight inetria term J33

uJ11            = ureal('uJ11',0,'Range',[-1,1]);                                       % uncertainty inertia term J11
uJ12            = ureal('uJ12',0,'Range',[-1,1]);                                       % uncertainty inertia term J12
uJ13            = ureal('uJ13',0,'Range',[-1,1]);                                       % uncertainty inertia term J13
uJ22            = ureal('uJ22',0,'Range',[-1,1]);                                       % uncertainty inertia term J22
uJ23            = ureal('uJ23',0,'Range',[-1,1]);                                       % uncertainty inertia term J23
uJ33            = ureal('uJ33',0,'Range',[-1,1]);                                       % uncertainty inertia term J33


Jnom            = [94.76,0.5929,1.416;0.5929,125.2,-1.283;1.416,-1.283,65.95];          % Nominal Inertia matrix
J               = [Jnom(1,1)*(1+wd11*uJ11),Jnom(1,2)*(1+wd12*uJ12),Jnom(1,3)*(1+wd13*uJ13);...
                   Jnom(2,1)*(1+wd12*uJ12),Jnom(2,2)*(1+wd22*uJ22),Jnom(2,3)*(1+wd23*uJ23);...
                   Jnom(3,1)*(1+wd13*uJ13),Jnom(3,2)*(1+wd23*uJ23),Jnom(3,3)*(1+wd33*uJ33)];
Jinv            = J^-1;

G               = ss(eye(3)/s^2);
G               = ss(G.a,G.b*Jinv,G.c,G.d);

% Define Reaction wheel model
Jrwl            = 0.00106;
Grwl            = ss(1/(Jrwl*s));
T               = 11.2;
w               = 2*pi/T;
damp            = 1.1;
alp             = 0.062;
Kp              = Jrwl*w^2;
Kd              = Jrwl*2*w*damp*0;
Ki              = Jrwl*2*alp*w^3;
Krwl            = balreal(ss((Kp + Ki/s)));
CLrwl           = zpk(feedback(Krwl*Grwl,1));

wRLWsf          = 0.05;                                                                 % RWL scale factor 5%
uRWL1sf         = ureal('uRWL1sf',0,'Range',[-1,1]);                                    % RWL1 scale factor uncertainty
uRWL2sf         = ureal('uRWL2sf',0,'Range',[-1,1]);                                    % RWL2 scale factor uncertainty
uRWL3sf         = ureal('uRWL3sf',0,'Range',[-1,1]);                                    % RWL3 scale factor uncertainty
uRWL4sf         = ureal('uRWL4sf',0,'Range',[-1,1]);                                    % RWL4 scale factor uncertainty

G_RWL           = blkdiag(CLrwl*(1+wRLWsf*uRWL1sf),CLrwl*(1+wRLWsf*uRWL2sf),CLrwl*(1+wRLWsf*uRWL3sf),CLrwl*(1+wRLWsf*uRWL4sf));

RLWbof          = [-0.5     ,-0.5     ,-0.5     ,-0.5;     ...                          % RWL mounting matrix
                    0.707206,-0.707206,-0.707206,0.707206; ...
                   -0.5     ,-0.5     , 0.5     ,0.5     ];
wRLWma          = 0.35/180*pi';                                                         % RWL misalignment
uDCM            = [null(RLWbof(:,1)')*wRLWma*[ ureal('uRWL1maX',0,'Range',[-1,1]);  ... % RWL1 misalignment x
                                               ureal('uRWL1maY',0,'Range',[-1,1])], ... % RWL1 misalignment y
                   null(RLWbof(:,2)')*wRLWma*[ ureal('uRWL2maX',0,'Range',[-1,1]);  ... % RWL2 misalignment x
                                               ureal('uRWL2maY',0,'Range',[-1,1])], ... % RWL2 misalignment y
                   null(RLWbof(:,3)')*wRLWma*[ ureal('uRWL3maX',0,'Range',[-1,1]);  ... % RWL3 misalignment x
                                               ureal('uRWL3maY',0,'Range',[-1,1])], ... % RWL3 misalignment y
                   null(RLWbof(:,4)')*wRLWma*[ ureal('uRWL4maX',0,'Range',[-1,1]);  ... % RWL4 misalignment x
                                               ureal('uRWL4maY',0,'Range',[-1,1])]];    % RWL4 misalignment y

RWLmat          = (RLWbof+uDCM)*G_RWL*pinv(RLWbof);                                     % overall RWL model

% Define controller
aK              = [-0.000005267527321,1,0,0; ...
                    0,-1.961613581091514,-0.477404761763747,1.327200434606469; ...
                    0,0,-0.592861929400114,1; ...
                    0,0,-0.620617904029083,-0.592861929400114];
bK              = [0;0;0;2.061773593367647];
cK              = [0.001317075749618,-0.621100505112127,-0.158399796833673,0.440356477431332];
dK              = 0;
K               = ss(aK,bK,cK,dK)*diag(diag(Jnom));

% Define generalized plant
systemnames     = 'G RWLmat';
inputvar        = '[r{3};n{3};u{3}]';
outputvar       = '[r-G;r-n-G]';
input_to_G      = '[RWLmat]';
input_to_RWLmat = '[u]';
cleanupsysic    = 'yes';
olic            = sysic;
P               = lft(olic(:,[1:3,7:9]),K);

% Perforance sensitivity weight
W               = (s+0.05)^2/2/(s+1e-5)^2*eye(3);

% Plot step response
figure;step(P)

% Compute the worst case gain for the full set of uncertainties
PW              = W*P;

% number of samples
N               = 500;

% Perform the sensitivity analysis with Morris method
out_morris      = sensanalysis(P,meth,N);

% Exclude the parameters that do affect the performance
PWred1          = usubs(PW,'uJ22',0,'uJ33',0,'uRWL1sf',0,'uRWL2sf',0,'uRWL3sf',0,'uRWL4sf',0);

% Exclude the parameters that do NOT affect the performance 
PWred2          = usubs(PW,'uJ11',0,'uJ12',0,'uJ13',0,'uJ23',0,'uRWL1maX',0,'uRWL1maY',0,'uRWL2maX',0,'uRWL2maY',0,'uRWL3maX',0,'uRWL3maY',0,'uRWL4maX',0,'uRWL4maY',0);

% Compute the worst case gains for different configurations
wcg_nom         = norm(PW.NominalValue,inf);
[wcg_red1,~]    = wcgain(PWred1);
[wcg_red2,~]    = wcgain(PWred2);
[wcg_tot,~]     = wcgain(PW);

disp(['Worst-case gain nominal system:             ', num2str(wcg_nom)]);
disp(['Worst-case gain reduced uncertain system 1: ', num2str(num2str(wcg_red1.UpperBound))]);
disp(['Worst-case gain reduced uncertain system 2: ', num2str(num2str(wcg_red2.UpperBound))]);
disp(['Worst-case gain complete uncertain system:  ', num2str(num2str(wcg_tot.UpperBound))]);