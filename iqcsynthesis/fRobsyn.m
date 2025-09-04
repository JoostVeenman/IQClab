function [K,ga] = fRobsyn(Pl,Delta,nout,nin,options)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        26-03-2020
% 
% -------------------------------------------------------------------------
%
% Description: This function performs a general IQC based robust controller
%              synthesis.
%
% Syntax:      [K,gamma]         = fRobsyn(Pl,Delta,nout,nin)
%              [K,gamma]         = fRobsyn(Pl,Delta,nout,nin,options)
%
% Usage:       This function solves the robust controller synthesis problem
%              based on IQCs by means of a D/K-type iterative controller
%              synthesis / robustness analysis algorithm.
%
%              A successful synthesis will yield a robust controller, which
%              guarantees robust stability and L2-gain performance for all
%              modelled uncertainties in the provided uncertainty set.
%
%              As input one should provide:
%
%                # The weighted uncertain open-loop generalized plant Pl
%                  with realization
%
%                             [ A  | Bp,  Bw,  Bu  ]
%                       [q]   [----|---------------] [p]
%                  Pl = [z] = [ Cq | Dqp, Dqw, Dqu ] [w]
%                       [y]   [ Cz | Dzp, Dzw, Dzu ] [u]
%                             [ Cy | Dyp, Dyw ,Dyu ]
%
%                  and with:
%                    - uncertainty input channel:  p
%                    - uncertainty output channel: q
%                    - disturbance input:          w 
%                    - performance output:         z 
%                    - control input:              u
%                    - measurement output:         y
%                    - (A,Bu) stabilizable 
%                    - (A,Cy) detectable
%
%                # The uncertainty block Delta = {Delta1,...,DeltaN},
%                  whose entries are objects from the class iqcdelta and
%                  which have been associated with an IQC multiplier using
%                  "iqcassign". 
%
%                  The uncertainty block Delta should be provided as a
%                  cell (see the function "iqcanalysis" for further
%                  information).
%
%                  Note: Even though it is not required by the iqcanalysis
%                  function, this function requires that the uncertainties
%                  and the corresponding input/output channels are ordered
%                  as Delta = {Delta_1, ... , Delta_N}.
%
%                # The input and output channel data:
%
%                  nin = [p;w;u]  and  nout = [q;z;y]
%
%                  where:
%
%                    - np denotes the number of uncertainty inputs
%                    - nw denotes the number of disturbance inputs
%                    - nu denotes the number of control inputs
%
%                  and where
%
%                    - nq denotes the number of uncertainty outputs
%                    - nz denotes the number of performance outputs
%                    - ny denotes the number of measurements outputs
%
%                  Note: As mentioned above, it is important that the
%                  open-loop plant obeys the order of the in- and output
%                  channels as suggested above.
%
%                # The optional structure "options", where:
%
%                   - options.maxiter     Maximum number of synthesis -
%                                         analysis iterations (default =
%                                         10).
%
%                   - options.subopt      If larger than 1, the algorithm
%                                         computes a suboptimal solution
%                                         options.subopt*gamma while
%                                         minimizing the norms on the
%                                         LMI variables (Default = 1.05).
%
%                   - options.condnr      If larger than 1, the algorithm
%                                         fixes the suboptimal solution
%                                         options.condnr*gamma while
%                                         minimizing the
%                                         conditioningsnumber of the
%                                         coupling condition
%                                         [Y22,I;I,X22] > 0 (default = 1)
%
%                   - options.constants   options.constants = [c1,c2,c3] 
%                                         is a verctor whose elements
%                                         perturb the LMIs with LMIi < -ciI
%                                         (Thus ci should be small)
%                                         (default = 1e-9*[1,1,1])
%
%                                         See the IQC-toolbox user-manual
%                                         for details on which LMI is
%                                         perturbed by which constant.
%
%                   - options.gmax        Specifies the maximal values of
%                                         gamma Might speed up the algoritm
%                                         (Default = 1000)
%
%                   - options.Pi11pos     Turns on the positivity condition
%
%                                         Pi11(iw) > options.Pi11pos*I
%
%                                         for all w\in R\cup{\infty}.
%
%                                         The default value is 1e-6. The
%                                         constraint is active if (small)
%                                         values larger than 0 are
%                                         specified.
%
%                                         Note, Pi11 is positive
%                                         semi-definite by nature. The
%                                         reason to add this constraint is
%                                         to enforce strict possitivity of
%                                         Pi11 which is needed to improve
%                                         the conditioning of the control
%                                         synthesis algorithm.
%
%                   - options.Parser      Use LMIlab or Yalmip to solve the
%                                         optimization problem. 
%
%                   - options.Solver      Use mincx in case LMIlab is
%                                         considered or another solver
%                                         (e.g. sedumi, sdpt3, mosek, etc.)
%                                         in case Yalmip is used.
%
%                   - options.FeasbRad    Bound variables (default = 1e9)
%
%                   - options.Terminate   can be used to change the
%                                         LMIlab-solver details (see help:
%                                         minxc) (default = 20)
%
%                   - options.RelAcc      can be used to change de
%                                         LMIlab-solver details (see help:
%                                         minxc) (default = 1e-4)
%
%              As Output, if feasible:
%
%                # For each synthesis-analysis iteration, the algorithm
%                  provides the robust controller u = K{i}y =
%                  ss(Ak{i},Bk{i},Ck{i},Dk{i}).
%
%                  For all Delta in the Delta set, these controllers
%                  stabilize the closed-loop system
%                  lft(lft(Delta,Pl),K{i}), while the L2-gain is rendered
%                  less than "ga(i)".
%
%                # The worst case L2-gain "ga(i)".
% -------------------------------------------------------------------------
if nargin == 4
    options.maxiter           = 10;
    options.subopt            = 1.01;
    options.condnr            = 1;
    options.Lyapext           = 3;
    options.constants         = 1e-9*ones(1,3);
    options.gmax              = 1000;    
    options.Pi11pos           = 1e-6;
    options.Parser            = 'LMIlab';
    options.Solver            = 'mincx';
    options.FeasbRad          = 1e9;
    options.Terminate         = 20;
    options.RelAcc            = 1e-4;
elseif nargin == 5
    if isfield(options,'maxiter')
        if options.maxiter >= 1
            options.maxiter   = options.maxiter;
        else
            options.maxiter   = 10;
        end
    else
        options.maxiter       = 10;
    end
    if isfield(options,'subopt')
        if options.subopt >= 1
            options.subopt    = options.subopt;
        else
            options.subopt    = 1.01;
        end
    else
        options.subopt        = 1.01;
    end
    if isfield(options,'condnr')
        if options.condnr >= 1
            options.condnr    = options.condnr;
        else
            options.condnr    = 1;
        end
    else
        options.condnr        = 1;
    end    
    if isfield(options,'Lyapext')
        if options.Lyapext == 1 || options.Lyapext == 2 || options.Lyapext == 3
            options.Lyapext   = options.Lyapext;
        else
            options.Lyapext   = 3;
        end
    else
        options.Lyapext       = 3;
    end
    if isfield(options,'constants')
        if min(options.constants) < 0
            options.constants = 1e-9*ones(1,3);
        else
            options.constants = options.constants;
        end
    else
        options.constants     = 1e-9*ones(1,3);
    end
    if isfield(options,'Parser')
        if strcmp(options.Parser,'LMIlab') || strcmp(options.Parser,'Yalmip')
            options.Parser    = options.Parser;
        else
            options.Parser    = 'LMIlab';
        end
    else
        options.Parser        = 'LMIlab';
    end
    if isfield(options,'Solver')
        if strcmp(options.Parser,'LMIlab')
            options.Solver    = 'mincx';
        elseif strcmp(options.Parser,'Yalmip')
            options.Solver    = options.Solver;
        end
    else
        options.Parser        = 'LMIlab';
        options.Solver        = 'mincx';
    end
    if isfield(options,'gmax')
        if options.gmax > 0
            options.gmax      = options.gmax;
        else
            options.gmax      = 1000;
        end
    else
        options.gmax          = 1000;
    end
    if isfield(options,'FeasbRad')
        if options.FeasbRad > 0
            options.FeasbRad  = options.FeasbRad;
        else
            options.FeasbRad  = 1e9;
        end
    else
        options.FeasbRad      = 1e9;
    end
    if isfield(options,'Terminate')
        if options.Terminate > 0
            options.Terminate = options.Terminate;
        else
            options.Terminate = 20;
        end
    else
        options.Terminate     = 20;
    end
    if isfield(options,'RelAcc')
        if options.RelAcc>0
            options.RelAcc    = options.RelAcc;
        else
            options.RelAcc    = 1e-4;
        end
    else
        options.RelAcc        = 1e-4;
    end
    if isfield(options,'Pi11pos')
        if options.Pi11pos > 0
            options.Pi11pos   = options.Pi11pos;
        else
            options.Pi11pos   = 1e-6;
        end
    else
        options.Pi11pos       = 1e-6;
    end
end

% Check input vector
if abs(size(Pl.d,2) - sum(nin)) > 0 || abs(size(Pl.d,1) - sum(nout)) > 0
    error('The the vectors "nin" or "nout" are not well-defined');
end

% Define dimensional info
np                            = nin(1);
nw                            = nin(2);
nu                            = nin(3);
nq                            = nout(1);
nz                            = nout(2);
ny                            = nout(3);

% Perform a nominal controller synthesis
clc
disp('-------------------------------------------------------------------');
disp(' ');
disp('Perform a nominal H-infinity synthesis...');
disp(' ');
disp('-------------------------------------------------------------------');

[K{1},ga_nom]                 = fWhinfsyn(Pl(nq+1:end,np+1:end),ss([]),[0,nz,ny],[0,nw,nu],options);

if ga_nom == -1
    error('LMIs are not feasible.');
else
    disp(['The nominal H-infinity synthesis problem is feasible with: gamma = ',num2str(ga_nom)])
end

% Perform a robustness analysis
disp('-------------------------------------------------------------------');
disp(' ');
disp('Perform a first IQC analysis...');
disp('-------------------------------------------------------------------');

N{1}                          = lft(Pl,K{1});
perf                          = iqcdelta('perf','ChannelClass','P','InputChannel',np+1:np+nw,'OutputChannel',nq+1:nq+nz,'PerfMetric','L2');
iqcOpt.FeasbRad               = 0.1*options.FeasbRad;
iqcOpt.Parser                 = options.Parser;
iqcOpt.Solver                 = options.Solver;
iqcOpt.Terminate              = options.Terminate;
iqcOpt.RelAcc                 = options.RelAcc;
iqcOpt.Pi11pos                = options.Pi11pos;
prob                          = iqcanalysis(N{1},{Delta,perf},iqcOpt);
ga(1)                         = prob.gamma;

if ga(1) == -1
    error('LMIs are not feasible.');
else
    disp('');
    disp(['The IQC analysis problem is feasible with: gamma = ',num2str(ga(1))])
end

% Perform the synthesis-analysis iteration
ga_new                        = ga(1);
ga_old                        = ga_nom;
i                             = 2;

while abs(ga_new-ga_old) > 1e-2 && i < options.maxiter
    % Perform a weighted controller synthesis
    P                         = iqcdec2mat(prob,prob.P);
    sc                        = prob.sc;
    Psi                       = prob.Psi;
    W                         = ssbal(Psi'*sc'*P*sc*Psi); % Construct the multiplier

    disp(' ');
    disp('Perform a weighted controller synthesis...');
    options.FeasbRad          = 10*options.FeasbRad;
    [K{i},ga(i)]              = fWhinfsyn(Pl,W,nout,nin,options);
    options.FeasbRad          = 0.1*options.FeasbRad;

    if ga(i) == -1
        error('LMIs are not feasible.');
    else
        disp(['The weighted synthesis problem is feasible with: gamma = ',num2str(ga(i))])
    end
    
    % Perform a robustness analysis
    disp(' ');
    disp('Perform a first IQC analysis...');
    N{i}                      = lft(Pl,K{i});
    prob                      = iqcanalysis(N{i},{Delta,perf},iqcOpt);
    ga(i)                     = prob.gamma;
    
    if ga(i) == -1
        error('LMIs are not feasible.');
    else
        disp(['The IQC analysis problem is feasible with: gamma = ',num2str(ga(i))])
    end
    ga_old                    = ga(i);
    ga_new                    = ga(i-1);
    i                         = i + 1;
end
end