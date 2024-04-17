function [K,ga] = fLPVsyn(Pl,Delta,nout,nin,options)
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
% Date:        22-03-2012
% 
% -------------------------------------------------------------------------
%
% Description: This function performs an LPV synthesis
%
% Syntax:      [K,gamma]         = fLPVsyn(Pl,Delta,nout,nin)
%              [K,gamma]         = fLPVsyn(Pl,Delta,nout,nin,options)
%
% Usage:       This function solves a particular version of the LPV
%              synthesis problem
%
%              As input:
%
%                # The plant Pl
%
%                              [ A  | Bs,  Bp,  Bu  ]
%                       [zs]   [----|---------------] [ws]
%                  Pl = [zp] = [ Cs | Dss, Dsp, Dsu ] [wz]
%                       [y ]   [ Cp | Dps, Dpp, Dpu ] [u ]
%                              [ Cy | Dys, Dyp ,Dyu ]
%
%                  with 
%                    - scheduling input:    ws
%                    - disturbance input:   wp
%                    - scheduling output:   zs 
%                    - performance output:  zp 
%                    - control input:       u
%                    - measurement output:  y
%                    - (A,Bu) stabilizable 
%                    - (A,Cy) detectable
%
%                # The uncertainty block Delta, which should be defined as
%                  an iqcdelta object, which must be compatible with the
%                  uncertainty class ultv: Delta is define as a linear map:
%
%
%                  Delta(delta) = sum_{i=1}^N delta_i*Ti = ...
%                                         = delta1*T1 + ... + deltaN*TN,
%
%                  where
%                    - Ti are some fixed matrices Ti \in R^{n x m}
%                    - delta:[0,\infty)->La is a piecewise continuous
%                      time-varying parameter vector that takes its values from
%                      the compact polytope
%
%                        La = co{delta^1, ... ,delta^M} = ...
%                                     = {sum_{a=1}^M b_a*delta^a: ...
%                                            b_a\leq0, sum_{a=1}^M b_a = 1}
%                                 
%                      with delta^j = (delta_1^j, ... ,delta_N^j),
%                      j\in{1,...,M}, as generator point and with 0\in La.
%                    - La is assumed to be star convex: [0,1]La\subset La.
%
%                # The input and output channel data:
%
%                  nin = [nws;nwp;nu]  and  nout = [nzs;nzp;ny]
%
%                  where:
%
%                    - nws denotes the number of scheduling inputs
%                    - nwp denotes the number of disturbance inputs
%                    - nu denotes the number of control inputs
%
%                  and where
%
%                    - nzs denotes the number of scheduling outputs
%                    - nzp denotes the number of performance outputs
%                    - ny denotes the number of measurements outputs
%
%                # The optional structure "options", where:
%
%                   - options.subopt      If larger than 1, the algorithm
%                                         computes a suboptimal solution
%                                         options.subopt*gamma while
%                                         minimizing the norms on the
%                                         LMI variables (default = 1.05).
%
%                   - options.condnr      If larger than 1, the algorithm
%                                         fixes the suboptimal solution
%                                         options.condnr*gamma while
%                                         minimizing the
%                                         conditioningsnumber of the
%                                         coupling condition [Y,I;I,X] > 0
%                                         (default = 1).
%
%                   - options.Lyapext     Constructs the extended Lyapunov
%                                         matrix "Xe" in different fashions
%                                         (default = 5, see fLyapext for
%                                         further details).
%
%                   - options.Multiplext  Constructs the extended
%                                         multiplier matrix "P_e" in
%                                         different fashions (Default = 3,
%                                         see fMultiplext for further
%                                         details).
%
%                   - options.RelaxType   Specifies the multiplier
%                                         relaxation type. Options are:
%                                           1.) 'DG': DG-scalings
%                                           2.) 'CH': Convex Hull (default)
%
%                   - options.constants   options.constants =
%                                                     = [c1,c2,c3,c4,c5,c6] 
%                                         is a verctor whose elements
%                                         perturb the LMIs with LMIi < -ciI
%                                         (Thus ci should be small)
%                                         (default = 1e-9*[1,1,1,1,1,1])
%
%                                         See the IQC-toolbox user-manual
%                                         for details on which LMI is
%                                         perturbed by which constant.
%
%                   - options.bounds      options.bounds = [b1,b2,b3,b4] 
%                                         is a verctor whose elements bound
%                                         the LMI variables if
%                                         options.subopt is set to values
%                                         larger than 1:
%
%                                         ||X||  < b1, ||Y||  < b2
%                                         ||Pp|| < b3 ,||Pd|| < b4
%
%                                         (default = [1,1,1,1])
%
%                   - options.Parser      Use LMIlab or Yalmip to solve the
%                                         optimization problem. 
%
%                   - options.Solver      Use mincx in case LMIlab is
%                                         considered or another solver
%                                         (e.g. sedumi, sdpt3, mosek, etc.)
%                                         in case Yalmip is used.
%
%                   - options.gmax        Specifies the maximal values of
%                                         gamma Might speed up the algoritm
%                                         (default = 1000)
%
%                   - options.FeasbRad    Bound variables (default = 1e9)
%
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
%                # The Controller u = lft(K,Delta_c(Delta))y, where K
%                  and Delta_c(Delta) are respectively defined as
%
%                              [ A | By,  Bw ]                              
%                  K =  [u ] = [---|---------] [y ]
%                       [zc]   [Cu | Duy, Duw] [wc]
%                              [Cc | Dcy, Dcw]
%
%                  and
%
%                  wc = Delta_c(Delta)zc = [0    , Delta^T] zc
%                                          [Delta,0       ]
%
%                  For all Delta\in bf{Delta}, this controller stabilizes
%                  the closed-loop system lft(Pl,lft(K,Delta_c(Delta))),
%                  while the following performance criterion is satisfied:
%
%                  int_0^\infty [zp(t)]^T[ I,       0 ][zp(t)]dt\geq0
%                               [zp(t)]  [ 0, -ga^2 I ][wp(t)]
%
%                  for all t\geq0.
%
%                # The L2-gain "ga".
% -------------------------------------------------------------------------
if nargin == 4
    options.subopt                 = 1.05;
    options.condnr                 = 1;
    options.Lyapext                = 5;
    options.Multiplext             = 3;
    options.RelaxType              = 'CH';
    options.constants              = 1e-9*ones(1,6);
    options.bounds                 = ones(1,4);
    options.Parser                 = 'LMIlab';
    options.Solver                 = 'mincx';
    options.gmax                   = 1000;
    options.FeasbRad               = 1e9;
    options.Terminate              = 20;
    options.RelAcc                 = 1e-4;
elseif nargin == 5
    if isfield(options,'subopt')
        if options.subopt >= 1
            options.subopt         = options.subopt;
        else
            options.subopt         = 1.05;
        end
    else
        options.subopt             = 1.05;
    end
    if isfield(options,'condnr')
        if options.condnr >= 1
            options.condnr         = options.condnr;
        else
            options.condnr         = 1;
        end
    else
        options.condnr             = 1;
    end    
    if isfield(options,'Lyapext')
        if options.Lyapext == 1 || options.Lyapext == 2 || options.Lyapext == 3 || options.Lyapext == 4 || options.Lyapext == 5
            options.Lyapext        = options.Lyapext;
        else
            options.Lyapext        = 5;
        end
    else
        options.Lyapext            = 5;
    end
    if isfield(options,'Multiplext')
        if options.Multiplext == 1 || options.Multiplext == 2 || options.Multiplext == 3
            options.Multiplext     = options.Multiplext;
        else
            options.Multiplext     = 3;
        end
    else
        options.Multiplext         = 3;
    end
    if isfield(options,'RelaxType')
        if strcmp(options.RelaxType,'DG') || strcmp(options.RelaxType,'CH') %|| strcmp(options.RelaxType,'PC') || strcmp(options.RelaxType,'ZP')
            options.RelaxType     = options.RelaxType;
        else
            options.RelaxType     = 'CH';
        end
    else
        options.RelaxType         = 'CH';
    end
    if isfield(options,'constants')
        if min(options.constants) < 0
            options.constants      = 1e-9*ones(1,6);
        else
            options.constants      = options.constants;
        end
    else
        options.constants          = 1e-9*ones(1,6);
    end
    if isfield(options,'bounds')
        if min(options.bounds) < 0
            options.bounds         = ones(1,4);
        else
            options.bounds         = options.bounds;
        end
    else
        options.bounds             = ones(1,4);
    end
    if isfield(options,'Parser')
        if strcmp(options.Parser,'LMIlab') || strcmp(options.Parser,'Yalmip')
            options.Parser         = options.Parser;
        else
            options.Parser         = 'LMIlab';
        end
    else
        options.Parser             = 'LMIlab';
    end
    if isfield(options,'Solver')
        if strcmp(options.Parser,'LMIlab')
            options.Solver         = 'mincx';
        elseif strcmp(options.Parser,'Yalmip')
            options.Solver         = options.Solver;
        end
    else
        options.Parser             = 'LMIlab';
        options.Solver             = 'mincx';
    end
    if isfield(options,'gmax')
        if options.gmax > 0
            options.gmax           = options.gmax;
        else
            options.gmax           = 1000;
        end
    else
        options.gmax               = 1000;
    end
    if isfield(options,'FeasbRad')
        if options.FeasbRad > 0
            options.FeasbRad       = options.FeasbRad;
        else
            options.FeasbRad       = 1e9;
        end
    else
        options.FeasbRad           = 1e9;
    end
    if isfield(options,'Terminate')
        if options.Terminate>0
            options.Terminate      = options.Terminate;
        else
            options.Terminate      = 20;
        end
    else
        options.Terminate          = 20;
    end
    if isfield(options,'RelAcc')
        if options.RelAcc>0
            options.RelAcc         = options.RelAcc;
        else
            options.RelAcc         = 1e-4;
        end
    else
        options.RelAcc             = 1e-4;
    end
end

% Check input vector
if abs(size(Pl.d,2) - sum(nin)) > 0 || abs(size(Pl.d,1) - sum(nout)) > 0
    error('The the vectors "nin" or "nout" are not well-defined');
end

% Define dimensional info
nws                                = nin(1);
nwp                                = nin(2);
nu                                 = nin(3);
nzs                                = nout(1);
nzp                                = nout(2);
ny                                 = nout(3);
nx                                 = size(Pl.a,1);

% Define ss matrices
A                                  = Pl.a;
Bs                                 = Pl.b(:,1:nws);
Bp                                 = Pl.b(:,nws+1:nws+nwp);
Bu                                 = Pl.b(:,nws+nwp+1:end);
Cs                                 = Pl.c(1:nzs,:);
Cp                                 = Pl.c(nzs+1:nzs+nzp,:);
Cy                                 = Pl.c(nzs+nzp+1:end,:);
Dss                                = Pl.d(1:nzs,1:nws);
Dsp                                = Pl.d(1:nzs,nws+1:nws+nwp);
Dsu                                = Pl.d(1:nzs,nws+nwp+1:end);
Dps                                = Pl.d(nzs+1:nzs+nzp,1:nws);
Dpp                                = Pl.d(nzs+1:nzs+nzp,nws+1:nws+nwp);
Dpu                                = Pl.d(nzs+1:nzs+nzp,nws+nwp+1:end);
Dys                                = Pl.d(nzs+nzp+1:end,1:nws);
Dyp                                = Pl.d(nzs+nzp+1:end,nws+1:nws+nwp);
Dyu                                = Pl.d(nzs+nzp+1:end,nws+nwp+1:end);

% create IQC problem
probOpt.FeasbRad                   = options.FeasbRad;
probOpt.Parser                     = options.Parser;
probOpt.Solver                     = options.Solver;
probOpt.Terminate                  = options.Terminate;
probOpt.RelAcc                     = options.RelAcc;
prob                               = iqcprob(probOpt);

if ~isobject(Delta)
    if ~strcmp(class(Delta),'iqcdelta')
        error('The scheduling parameter block Delta should be defined as an object from the uncertainty class iqcdelta')
    end
else
    % Create primal full-block multiplier 
    prob.eps                       = options.constants(5);
    uDelta_p                       = iqcassign(Delta,'ultv','RelaxationType',options.RelaxType,'PrimalDual','Primal');
    prob                           = iqcltv_p(uDelta_p,prob);
    
    % Create dual full-block multiplier
    prob.eps                       = options.constants(6);
    uDelta_d                       = iqcassign(Delta,'ultv','RelaxationType',options.RelaxType,'PrimalDual','Dual');
    prob                           = iqcltv_d(uDelta_d,prob);
end

prob.eps                           = 0;

Pp                                 = prob.P(1:nws+nzs,1:nws+nzs);
Pd                                 = prob.P(nws+nzs+1:2*(nws+nzs),nws+nzs+1:2*(nws+nzs));

% Define LMI variables
X                                  = iqcvar(prob,[nx,nx],'symmetric');
Y                                  = iqcvar(prob,[nx,nx],'symmetric');
gamma                              = iqcvar(prob,[1,1],'symmetric');

% Construct matrix variables
Vx                                 = blkdiag(oblkdiag(X),Pp,-gamma*eye(nzp+nwp));
Vy                                 = blkdiag(oblkdiag(Y),Pd, gamma*eye(nwp+nzp));
Vc                                 = blkdiag(Y,X);

% Construct outer factors
Top                                = [fJ(nx,nws+nwp+nzp)';A,Bs,Bp,zeros(nx,nzp);Cs,Dss,Dsp,zeros(nzs,nzp);fJ(nws,nwp+nzp,nx)';fJt(nwp+nzp,nx+nws)'];
Tod                                = [-A',-Cs',-Cp',zeros(nx,nwp);fJ(nx,nzs+nzp+nwp)';fJ(nzs,nwp+nzp,nx)';-Bs',-Dss',-Dps',zeros(nws,nwp);fJt(nzp+nwp,nx+nzs)'];
Tcp                                = -fHe(fJt(nzp,nx+nws+nwp)*[Cp,Dps,Dpp,zeros(nzp)]);
Tcd                                = -fHe(fJt(nwp,nx+nzs+nzp)*[-Bp',-Dsp',-Dpp',zeros(nwp)]);

% Construct annihilators
Tap                                = blkdiag(null([Cy,Dys,Dyp]),eye(nzp));
Tad                                = blkdiag(null([Bu',Dsu',Dpu']),eye(nwp));

% Perturb LMIs by small delta
Ep                                 = -options.constants(3)*eye(size(Tap'*Tap,1));
Ed                                 = options.constants(4)*eye(size(Tad'*Tad,1));
Ec                                 = blkdiag((options.constants(1))*eye(nx),(options.constants(2))*eye(nx));
    
% Construct LMI
prob                               = iqclmi(prob,Vx,-1,Tap'*Tcp*Tap+Ep,Top*Tap);
prob                               = iqclmi(prob,Vy, 1,Tad'*Tcd*Tad+Ed,Tod*Tad);
prob                               = iqclmi(prob,Vc, 1,-fOblkdiag(eye(nx))+Ec);
prob                               = iqclmi(prob,gamma,-1,options.gmax);

% solve LPVsyn problem
prob                               = iqcsolve(prob,gamma);
ga                                 = prob.gamma;

if ga == -1
    error('LMIs are not feasible.');
else
    disp(['The LPV synthesis problem is feasible with: gamma = ',num2str(ga)])
end

X                                  = iqcdec2mat(prob,X);
Y                                  = iqcdec2mat(prob,Y);
Pp                                 = iqcdec2mat(prob,Pp);
Pd                                 = iqcdec2mat(prob,Pd);

disp(' ');
disp(['The norm of X  is ', num2str(norm(X),'%10.4e\n')]);
disp(['The norm of Y  is ', num2str(norm(Y),'%10.4e\n')]);
disp(['The norm of Pp is ', num2str(norm(Pp),'%10.4e\n')]);
disp(['The norm of Pd is ', num2str(norm(Pd),'%10.4e\n')]);

% Compute suboptimal solution, while minimizing the norm of the variables
if options.subopt > 1
    % create IQC problem
    probOpt.FeasbRad               = options.FeasbRad;
    probOpt.Parser                 = options.Parser;
    probOpt.Solver                 = options.Solver;
    probOpt.Terminate              = options.Terminate;
    probOpt.RelAcc                 = options.RelAcc;
    prob                           = iqcprob(probOpt);
    
    if ~isobject(Delta)
        if ~strcmp(class(Delta),'iqcdelta')
            error('The scheduling parameter block Delta should be defined as an object from the uncertainty class iqcdelta')
        end
    else
        % Create primal full-block multiplier 
        prob.eps                   = options.constants(5);
        uDelta_p                   = iqcassign(Delta,'ultv','RelaxationType',options.RelaxType,'PrimalDual','Primal');
        prob                       = iqcltv_p(uDelta_p,prob);

        % Create dual full-block multiplier
        prob.eps                   = options.constants(6);
        uDelta_d                   = iqcassign(Delta,'ultv','RelaxationType',options.RelaxType,'PrimalDual','Dual');
        prob                       = iqcltv_d(uDelta_d,prob);
    end

    prob.eps                       = 0;

    Pp                             = prob.P(1:nws+nzs,1:nws+nzs);
    Pd                             = prob.P(nws+nzs+1:2*(nws+nzs),nws+nzs+1:2*(nws+nzs));

    % Define LMI variables
    X                              = iqcvar(prob,[nx,nx],'symmetric');
    Y                              = iqcvar(prob,[nx,nx],'symmetric');
    gamma                          = iqcvar(prob,[1,1],'symmetric');

    % Bound LMIvariables
    prob                           = fBoundvar(prob,X,gamma,options.bounds(1));
    prob                           = fBoundvar(prob,Y,gamma,options.bounds(2));
    prob                           = fBoundvar(prob,Pp,gamma,options.bounds(3));
    prob                           = fBoundvar(prob,Pd,gamma,options.bounds(4));

    % Construct matrix variables
    Vx                             = blkdiag(oblkdiag(X),Pp);
    Vy                             = blkdiag(oblkdiag(Y),Pd);
    Vc                             = blkdiag(Y,X);

    % Construct outer factors
    ga                             = options.subopt*ga;
    Top                            = [fJ(nx,nws+nwp+nzp)';A,Bs,Bp,zeros(nx,nzp);Cs,Dss,Dsp,zeros(nzs,nzp);fJ(nws,nwp+nzp,nx)'];
    Tod                            = [-A',-Cs',-Cp',zeros(nx,nwp);fJ(nx,nzs+nzp+nwp)';fJ(nzs,nwp+nzp,nx)';-Bs',-Dss',-Dps',zeros(nws,nwp)];    
    Tcp                            = -fHe(fJt(nzp,nx+nws+nwp)*[Cp,Dps,Dpp,zeros(nzp)])+ga*fJt(nwp+nzp,nx+nws)*fJt(nwp+nzp,nx+nws)';
    Tcd                            = -fHe(fJt(nwp,nx+nzs+nzp)*[-Bp',-Dsp',-Dpp',zeros(nwp)])-ga*fJt(nzp+nwp,nx+nzs)*fJt(nzp+nwp,nx+nzs)';

    % Construct annihilators
    Tap                            = blkdiag(null([Cy,Dys,Dyp]),eye(nzp));
    Tad                            = blkdiag(null([Bu',Dsu',Dpu']),eye(nwp));

    % Perturb LMIs by small delta
    Ep                             = -options.constants(3)*eye(size(Tap'*Tap,1));
    Ed                             = options.constants(4)*eye(size(Tad'*Tad,1));
    Ec                             = blkdiag((options.constants(1))*eye(nx),(options.constants(2))*eye(nx));
    
    % Construct LMI
    prob                           = iqclmi(prob,Vx,-1,Tap'*Tcp*Tap+Ep,Top*Tap);
    prob                           = iqclmi(prob,Vy, 1,Tad'*Tcd*Tad+Ed,Tod*Tad);
    prob                           = iqclmi(prob,Vc, 1,-fOblkdiag(eye(nx))+Ec);
    prob                           = iqclmi(prob,gamma,1);
    
    % solve LPVsyn problem
    disp(' ');
    disp('Minimize the norm of the LMI variables...');
    prob                           = iqcsolve(prob,gamma);
    gambnd                         = prob.gamma;

    if gambnd == -1
         error('Suboptimal (2nd) LMIs are not feasible!');
    else
        disp(['The subobtial LPV synthesis problem is feasible with gamma = ',num2str(ga)])
    end
    
    X                              = iqcdec2mat(prob,X);
    Y                              = iqcdec2mat(prob,Y);
    Pp                             = iqcdec2mat(prob,Pp);
    Pd                             = iqcdec2mat(prob,Pd);

    disp(' ');
    disp(['The norm of X  is ', num2str(norm(X),'%10.4e\n')]);
    disp(['The norm of Y  is ', num2str(norm(Y),'%10.4e\n')]);
    disp(['The norm of Pp is ', num2str(norm(Pp),'%10.4e\n')]);
    disp(['The norm of Pd is ', num2str(norm(Pd),'%10.4e\n')]);
end

% Improve the conditioning number of the coupling condition
if options.condnr > 1
    if options.subopt > 1
        gambnd_subopt              = options.condnr*gambnd;
    else
        options.FeasbRad           = options.condnr*options.FeasbRad;
        gambnd_subopt              = options.condnr*options.FeasbRad;
        ga                         = options.condnr*ga;
    end

    % create IQC problem
    probOpt.FeasbRad               = options.FeasbRad;
    probOpt.Parser                 = options.Parser;
    probOpt.Solver                 = options.Solver;
    probOpt.Terminate              = options.Terminate;
    probOpt.RelAcc                 = options.RelAcc;
    prob                           = iqcprob(probOpt);
    
    if ~isobject(Delta)
        if ~strcmp(class(Delta),'iqcdelta')
            error('The scheduling parameter block Delta should be defined as an object from the uncertainty class iqcdelta')
        end
    else
        % Create primal full-block multiplier 
        prob.eps                   = options.constants(5);
        uDelta_p                   = iqcassign(Delta,'ultv','RelaxationType',options.RelaxType,'PrimalDual','Primal');
        prob                       = iqcltv_p(uDelta_p,prob);

        % Create dual full-block multiplier
        prob.eps                   = options.constants(6);
        uDelta_d                   = iqcassign(Delta,'ultv','RelaxationType',options.RelaxType,'PrimalDual','Dual');
        prob                       = iqcltv_d(uDelta_d,prob);
    end

    prob.eps                       = 0;

    Pp                             = prob.P(1:nws+nzs,1:nws+nzs);
    Pd                             = prob.P(nws+nzs+1:2*(nws+nzs),nws+nzs+1:2*(nws+nzs));

    % Define LMI variables
    X                              = iqcvar(prob,[nx,nx],'symmetric');
    Y                              = iqcvar(prob,[nx,nx],'symmetric');
    gamma                          = iqcvar(prob,[1,1],'symmetric');

    % Bound LMIvariables
    prob                           = fBoundvar(prob,X,gambnd_subopt,options.bounds(1));
    prob                           = fBoundvar(prob,Y,gambnd_subopt,options.bounds(2));
    prob                           = fBoundvar(prob,Pp,gambnd_subopt,options.bounds(3));
    prob                           = fBoundvar(prob,Pd,gambnd_subopt,options.bounds(4));
    
    % Construct matrix variables
    Vx                             = blkdiag(oblkdiag(X),Pp);
    Vy                             = blkdiag(oblkdiag(Y),Pd);
    Vc                             = blkdiag(Y,X);
    Txy                            = [Y,-gamma*eye(nx);-gamma*eye(nx),X];

    % Construct outer factors
    Top                            = [fJ(nx,nws+nwp+nzp)';A,Bs,Bp,zeros(nx,nzp);Cs,Dss,Dsp,zeros(nzs,nzp);fJ(nws,nwp+nzp,nx)'];
    Tod                            = [-A',-Cs',-Cp',zeros(nx,nwp);fJ(nx,nzs+nzp+nwp)';fJ(nzs,nwp+nzp,nx)';-Bs',-Dss',-Dps',zeros(nws,nwp)]; 
    Tcp                            = -fHe(fJt(nzp,nx+nws+nwp)*[Cp,Dps,Dpp,zeros(nzp)])+ga*fJt(nwp+nzp,nx+nws)*fJt(nwp+nzp,nx+nws)';
    Tcd                            = -fHe(fJt(nwp,nx+nzs+nzp)*[-Bp',-Dsp',-Dpp',zeros(nwp)])-ga*fJt(nzp+nwp,nx+nzs)*fJt(nzp+nwp,nx+nzs)';
    
    % Construct annihilators
    Tap                            = blkdiag(null([Cy,Dys,Dyp]),eye(nzp));
    Tad                            = blkdiag(null([Bu',Dsu',Dpu']),eye(nwp));

    % Perturb LMIs by small delta
    Ep                             = -options.constants(3)*eye(size(Tap'*Tap,1));
    Ed                             = options.constants(4)*eye(size(Tad'*Tad,1));
    Ec                             = blkdiag((options.constants(1))*eye(nx),(options.constants(2))*eye(nx));
    
    % Construct LMI
    prob                           = iqclmi(prob,Vx,-1,Tap'*Tcp*Tap+Ep,Top*Tap);
    prob                           = iqclmi(prob,Vy, 1,Tad'*Tcd*Tad+Ed,Tod*Tad);
    prob                           = iqclmi(prob,Vc, 1,-fOblkdiag(eye(nx))+Ec);
    prob                           = iqclmi(prob,Txy,1,-fOblkdiag(eye(nx)));
    
    % solve LPVsyn problem
    disp(' ');
    disp('Improve the conditioning of the coupling condition...');
    prob                           = iqcsolve(prob,gamma);
    gam_cond                       = prob.gamma;

    if gam_cond == -1
         error('Suboptimal (3nd) LMIs are not feasible.');
    else
       disp(['The subobtial LPV synthesis problem is feasible with gamma = ',num2str(ga)]) 
    end

    % call variables
    X                              = iqcdec2mat(prob,X);
    Y                              = iqcdec2mat(prob,Y);
    Pp                             = iqcdec2mat(prob,Pp);
    Pd                             = iqcdec2mat(prob,Pd);
    
    disp(' ');
    disp(['The norm of X  is ', num2str(norm(X),'%10.4e\n')]);
    disp(['The norm of Y  is ', num2str(norm(Y),'%10.4e\n')]);
    disp(['The norm of Pp is ', num2str(norm(Pp),'%10.4e\n')]);
    disp(['The norm of Pd is ', num2str(norm(Pd),'%10.4e\n')]);
end
disp('--------------------------------------------------------------------')

% Construct the extended Lyapunov matrix
Xe                                 = fLyapext(X,Y,options.Lyapext);

% Construct the extended multiplier
[~,Qe,Se,Re]                       = fMultiplext(Pp,Pd,options.Multiplext);

% Construct the control function K
nzc                                = size(Qe,1)-nzs;
nwc                                = size(Re,1)-nws;

Qei                                = fInv(Qe);
Te                                 = blkdiag(0.5*fHe(Re-Se'*Qei*Se),-ga*eye(nwp));

J                                  = fJ(2*nx,nws+nwc+nwp)';
Jt                                 = fJt(nws+nwc+nwp,2*nx)';

A1                                 = Xe*[blkdiag(A,zeros(nx)),blkdiag(Bs,zeros(nx,nwc)),[Bp;zeros(nx,nwp)]];
B1                                 = Xe*[zeros(nx),Bu,zeros(nx,nzc);fJ(nx,nu+nzc)'];
A2                                 = [blkdiag(Cs,zeros(nzc,nx)),blkdiag(Dss,zeros(nzc,nwc))+Qe\Se,[Dsp;zeros(nzc,nwp)];Cp,zeros(nzp,nx),Dps,zeros(nzp,nwc),Dpp];
B2                                 = [zeros(nzs+nzc+nzp,nx),[blkdiag(Dsu,eye(nzc));[Dpu,zeros(nzp,nzc)]]];

T1                                 = [J'*B1;B2]';
T2                                 = [[fJ(nx,nws+nwc+nwp,nx)';Cy,zeros(ny,nx),Dys,zeros(ny,nwc),Dyp;fJ(nwc,nwp,2*nx+nws)'],zeros(nx+ny+nwc,nzs+nzc+nzp)];
T3                                 = [Jt'*Te*Jt+fHe(J'*A1),A2';A2,-blkdiag(Qei,ga*eye(nzp))];

Kt                                 = fInvproj(T1,T2,T3);
K                                  = (eye(size(Kt,1))+Kt*blkdiag(zeros(nx),Dyu,zeros(nzc,nwc)))\Kt;
Ak                                 = K(1:nx,1:nx);
Bk                                 = K(1:nx,nx+1:end);
Ck                                 = K(nx+1:end,1:nx);
Dk                                 = K(nx+1:end,nx+1:end);
K                                  = ss(Ak,Bk,Ck,Dk);

% Cosntruct the scheduling function 
V                                  = fInv(Re);
U                                  = Qe-Se*V*Se';
W                                  = V*Se';

U11                                = U(1:nzs,1:nzs);
U12                                = U(1:nzs,nzs+1:end);
U22                                = U(nzs+1:end,nzs+1:end);

V11                                = V(1:nws,1:nws);
V12                                = V(1:nws,nws+1:end);
V22                                = V(nws+1:end,nws+1:end);

W11                                = W(1:nws,1:nzs);
W12                                = W(1:nws,nzs+1:end);
W21                                = W(nws+1:end,1:nzs);
W22                                = W(nws+1:end,nzs+1:end);

A11                                = -W22;
A12                                = [W21,V12'];
A21                                = [U12;W12];
A22                                = [U11,W11';W11,V11];

fScheduling                        = ss([A11+A12/A22*A21,-A12/A22;A22\A21,-fInv(A22)]);
K                                  = lft(K,fScheduling,nwc,nzc);
end