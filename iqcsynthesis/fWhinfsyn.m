function [K,ga] = fWhinfsyn(Pl,W,nout,nin,options)
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
% Description: This function performs an weighted nominal controller
%              synthesis. This function is primarily used in fRobsyn.
%
% Syntax:      [K,gamma]         = fWhinfsyn(Pl,W,nout,nin)
%              [K,gamma]         = fWhinfsyn(Pl,W,nout,nin,options)
%
% Usage:       This function solves a particular version of the H-infinity
%              synthesis problem (allowing for a control synthesis with
%              unstable weights. 
%
%              Note: This function is used in the function fRobsyn to
%              perform the nominal (weighted) controller synthesis. See
%              the function fRobsyn for further details.
%
%              As input:
%
%                # The plant Pl
%
%                             [ A  | Bp,  Bw,  Bu  ]
%                       [q]   [----|---------------] [p]
%                  Pl = [z] = [ Cq | Dqp, Dqw, Dqu ] [w]
%                       [y]   [ Cz | Dzp, Dzw, Dzu ] [u]
%                             [ Cy | Dyp, Dyw ,Dyu ]
%
%                  with 
%                    - uncertainty input channel:  p
%                    - uncertainty output channel: q
%                    - disturbance input:          w 
%                    - performance output:         z 
%                    - control input:              u
%                    - measurement output:         y
%                    - (A,Bu) stabilizable 
%                    - (A,Cy) detectable
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
%                # The (unstable) weight W = ss(Aw,Bw,Cw,Dw)\in RL_\infty
%
%                  with
%                    - input/output dimension np + nq
%                    - no poles on the imaginary axis
%
%                # The optional structure "options", where:
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
%                                         (Default = 1000)
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
%                # The Controller u = Ky = ss(Ak,Bk,Ck,Dk).
%
%                  This controller stabilizes the closed-loop system
%                  lft(Pl,K), while the following performance criterion is
%                  satisfied: 
%
%                  int_0^\infty [z(t)]^T[ I,       0 ][z(t)]dt\geq0
%                               [w(t)]  [ 0, -ga^2 I ][w(t)]
%
%                  for all t\geq0.
%
%                # The L2-gain "ga".
% -------------------------------------------------------------------------
if nargin == 4
    options.subopt                 = 1.05;
    options.condnr                 = 1;
    options.Lyapext                = 3;
    options.constants              = 1e-9*ones(1,3);
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
        if options.Lyapext == 1 || options.Lyapext == 2 || options.Lyapext == 3
            options.Lyapext        = options.Lyapext;
        else
            options.Lyapext        = 3;
        end
    else
        options.Lyapext            = 3;
    end
    if isfield(options,'constants')
        if min(options.constants) < 0
            options.constants      = 1e-9*ones(1,3);
        else
            options.constants      = options.constants;
        end
    else
        options.constants          = 1e-9*ones(1,3);
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
np                                 = nin(1);
nw                                 = nin(2);
nu                                 = nin(3);
nq                                 = nout(1);
nz                                 = nout(2);
ny                                 = nout(3);
nx                                 = size(Pl.a,1);

% Define ss matrices
A                                  = Pl.a;
Bp                                 = Pl.b(:,1:np);
Bw                                 = Pl.b(:,np+1:np+nw);
Bu                                 = Pl.b(:,np+nw+1:end);

Cq                                 = Pl.c(1:nq,:);
Cz                                 = Pl.c(nq+1:nq+nz,:);
Cy                                 = Pl.c(nq+nz+1:end,:);

Dqp                                = Pl.d(1:nq,1:np);
Dqw                                = Pl.d(1:nq,np+1:np+nw);
Dqu                                = Pl.d(1:nq,np+nw+1:end);

Dzp                                = Pl.d(nq+1:nq+nz,1:np);
Dzw                                = Pl.d(nq+1:nq+nz,np+1:np+nw);
Dzu                                = Pl.d(nq+1:nq+nz,np+nw+1:end);

Dyp                                = Pl.d(nq+nz+1:end,1:np);
Dyw                                = Pl.d(nq+nz+1:end,np+1:np+nw);
Dyu                                = Pl.d(nq+nz+1:end,np+nw+1:end);

if size(W.b,2) ~= np + nq && size(W.c,1) ~= np + nq
   error('The dimension of the weighting function W should be equal to the number of uncertainty inputs plus the number of uncertainty outputs of Pl (i.e. np + nq).'); 
end

B                                  = [zeros(nx,nq),Bp,zeros(nx,nz),Bw];
B1                                 = Bu;
C                                  = 2*[Cq;zeros(np,nx);Cz];
C1                                 = Cy;
D                                  = 2*[-0.5*eye(nq),Dqp,zeros(nq,nz),Dqw;zeros(np,nq),0.5*eye(np),zeros(np,nz+nw);zeros(nz,nq),Dzp,zeros(nz,nz),Dzw];
E                                  = 2*[Dqu;zeros(np,nu);Dzu];
F                                  = [zeros(ny,nq),Dyp,zeros(ny,nz),Dyw];

nwx                                = size(W.a,1);
AW                                 = W.a;
BW                                 = [W.b,zeros(nwx,nz+nw)];
CW                                 = [W.c;zeros(nz,nwx)];
DW                                 = blkdiag(W.d,[eye(nz),zeros(nz,nw)]);
nwi                                = size(BW,2);
nwo                                = size(CW,1);

% create IQC problem
probOpt.FeasbRad                   = options.FeasbRad;
probOpt.Parser                     = options.Parser;
probOpt.Solver                     = options.Solver;
probOpt.Terminate                  = options.Terminate;
probOpt.RelAcc                     = options.RelAcc;
probOpt.eps                        = 0;
prob                               = iqcprob(probOpt);

% Define LMI variables
X                                  = iqcvar(prob,[nx+nwx,nx+nwx],'symmetric');
Y                                  = iqcvar(prob,[nx+nwx,nx+nwx],'symmetric');
gamma                              = iqcvar(prob,[1,1],'symmetric');

% Construct matrix variables
Mga                                = blkdiag(zeros(nq+np),-gamma*eye(nz+nw));
Vx                                 = blkdiag(oblkdiag(X),Mga);
Vy                                 = blkdiag(oblkdiag(Y),-Mga);
Vc                                 = blkdiag(Y,X);

% Construct outer factors
Top                                = [eye(nwx+nx),zeros(nwx+nx,nwi);blkdiag(AW,A),[BW;B];sqrt(2)*fJt(nwi,nwx+nx)'];
Tod                                = [eye(nwx+nx),zeros(nwx+nx,nwi);AW,zeros(nwx,nx),BW;-C'*CW,-A',-C'*DW;sqrt(2)*fJt(nwi,nwx+nx)'];
Toc                                = [fJ(nwx+nx,nx)';fJ(nwx,2*nx)';-fJt(nx,nwx+nx)'];
Tcp                                = -fHe([CW';zeros(nx,nwo);DW']*[zeros(nwo,nwx),C,D]);
Tcd                                = -fHe([-CW'*D;-B;-DW'*D]*fJt(nwi,nwx+nx)');
Tcc                                = blkdiag(zeros(nwx),fOblkdiag(-eye(nx)));

% Construct annihilators
Tap                                = null([zeros(size(C1,1),nwx),C1,F]);
Tad                                = null([E'*CW,B1',E'*DW]);

% Perturb LMIs by small delta
Ep                                 = -options.constants(2)*eye(size(Tap'*Tap,1));
Ed                                 = options.constants(3)*eye(size(Tad'*Tad,1));
Ec                                 = options.constants(1)*eye(nwx+2*nx);
    
% Construct LMI
prob                               = iqclmi(prob,Vx,-1,Tap'*Tcp*Tap+Ep,Top*Tap);
prob                               = iqclmi(prob,Vy, 1,Tad'*Tcd*Tad+Ed,Tod*Tad);
prob                               = iqclmi(prob,Vc, 1,Tcc+Ec,Toc);
prob                               = iqclmi(prob,gamma,-1,options.gmax);

% solve synthesis problem
prob                               = iqcsolve(prob,gamma);
ga                                 = prob.gamma;

if ga == -1
    error('LMIs are not feasible.');
else
    disp(['The synthesis problem is feasible with: gamma = ',num2str(ga)])
end

X                                  = iqcdec2mat(prob,X);
Y                                  = iqcdec2mat(prob,Y);

disp(' ');
disp(['The norm of X  is ', num2str(norm(X),'%10.4e\n')]);
disp(['The norm of Y  is ', num2str(norm(Y),'%10.4e\n')]);

% Compute suboptimal solution, while minimizing the norm of the variables
if options.subopt > 1
    % create IQC problem
    probOpt.FeasbRad               = options.FeasbRad;
    probOpt.Parser                 = options.Parser;
    probOpt.Solver                 = options.Solver;
    probOpt.Terminate              = options.Terminate;
    probOpt.RelAcc                 = options.RelAcc;
    probOpt.eps                    = 0;
    prob                           = iqcprob(probOpt);

    % Define LMI variables
    X                              = iqcvar(prob,[nx+nwx,nx+nwx],'symmetric');
    Y                              = iqcvar(prob,[nx+nwx,nx+nwx],'symmetric');
    gamma                          = iqcvar(prob,[1,1],'symmetric');

    % Bound LMIvariables
    prob                           = fBoundvar(prob,X,gamma);
    prob                           = fBoundvar(prob,Y,gamma);

    % Construct matrix variables
    Vx                             = oblkdiag(X);
    Vy                             = oblkdiag(Y);
    Vc                             = blkdiag(Y,X);    

    % Construct outer factors
    ga                             = options.subopt*ga;
    Mga                            = blkdiag(zeros(nq+np),-ga*eye(nz+nw));
    Top                            = [eye(nwx+nx),zeros(nwx+nx,nwi);blkdiag(AW,A),[BW;B]];
    Tod                            = [eye(nwx+nx),zeros(nwx+nx,nwi);AW,zeros(nwx,nx),BW;-C'*CW,-A',-C'*DW];
    Toc                            = [fJ(nwx+nx,nx)';fJ(nwx,2*nx)';-fJt(nx,nwx+nx)'];
    Tcp                            = -fHe([CW';zeros(nx,nwo);DW']*[zeros(nwo,nwx),C,D])-2*fJt(nwi,nwx+nx)*Mga*fJt(nwi,nwx+nx)';
    Tcd                            = -fHe([-CW'*D;-B;-DW'*D]*fJt(nwi,nwx+nx)')+2*fJt(nwi,nwx+nx)*Mga*fJt(nwi,nwx+nx)';
    Tcc                            = blkdiag(zeros(nwx),fOblkdiag(-eye(nx)));

    % Construct annihilators
    Tap                            = null([zeros(size(C1,1),nwx),C1,F]);
    Tad                            = null([E'*CW,B1',E'*DW]);

    % Perturb LMIs by small delta
    Ep                             = -options.constants(2)*eye(size(Tap'*Tap,1));
    Ed                             = options.constants(3)*eye(size(Tad'*Tad,1));
    Ec                             = options.constants(1)*eye(nwx+2*nx);
    
    % Construct LMI
    prob                           = iqclmi(prob,Vx,-1,Tap'*Tcp*Tap+Ep,Top*Tap);
    prob                           = iqclmi(prob,Vy, 1,Tad'*Tcd*Tad+Ed,Tod*Tad);
    prob                           = iqclmi(prob,Vc, 1,Tcc+Ec,Toc);
    
    % solve synthesis problem
    disp(' ');
    disp('Minimize the norm of the LMI variables...');
    prob                           = iqcsolve(prob,gamma);
    gambnd                         = prob.gamma;

    if gambnd == -1
         error('Suboptimal (2nd) LMIs are not feasible!');
    else
        disp(['The subobtial synthesis problem is feasible with gamma = ',num2str(ga)])
    end
    
    X                              = iqcdec2mat(prob,X);
    Y                              = iqcdec2mat(prob,Y);

    disp(' ');
    disp(['The norm of X  is ', num2str(norm(X),'%10.4e\n')]);
    disp(['The norm of Y  is ', num2str(norm(Y),'%10.4e\n')]);
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
    probOpt.eps                    = 0;
    prob                           = iqcprob(probOpt);

    % Define LMI variables
    X                              = iqcvar(prob,[nx+nwx,nx+nwx],'symmetric');
    Y                              = iqcvar(prob,[nx+nwx,nx+nwx],'symmetric');
    gamma                          = iqcvar(prob,[1,1],'symmetric');

    % Bound LMIvariables
    prob                           = fBoundvar(prob,X,gambnd_subopt);
    prob                           = fBoundvar(prob,Y,gambnd_subopt);
    
    % Construct matrix variables
    Vx                             = oblkdiag(X);
    Vy                             = oblkdiag(Y);
    Vc                             = blkdiag(Y,X);
    Vxy                            = [Y(nwx+1:nwx+nx,nwx+1:nwx+nx),-gamma*eye(nx);-gamma*eye(nx),X(nwx+1:nwx+nx,nwx+1:nwx+nx)];

    % Construct outer factors
    ga                             = options.condnr*ga;
    Mga                            = blkdiag(zeros(nq+np),-ga*eye(nz+nw));        
    Top                            = [eye(nwx+nx),zeros(nwx+nx,nwi);blkdiag(AW,A),[BW;B]];
    Tod                            = [eye(nwx+nx),zeros(nwx+nx,nwi);AW,zeros(nwx,nx),BW;-C'*CW,-A',-C'*DW];
    Toc                            = [fJ(nwx+nx,nx)';fJ(nwx,2*nx)';-fJt(nx,nwx+nx)'];
    Tcp                            = -fHe([CW';zeros(nx,nwo);DW']*[zeros(nwo,nwx),C,D])-2*fJt(nwi,nwx+nx)*Mga*fJt(nwi,nwx+nx)';
    Tcd                            = -fHe([-CW'*D;-B;-DW'*D]*fJt(nwi,nwx+nx)')+2*fJt(nwi,nwx+nx)*Mga*fJt(nwi,nwx+nx)';
    Tcc                            = blkdiag(zeros(nwx),fOblkdiag(-eye(nx)));

    % Construct annihilators
    Tap                            = null([zeros(size(C1,1),nwx),C1,F]);
    Tad                            = null([E'*CW,B1',E'*DW]);

    % Perturb LMIs by small delta
    Ep                             = -options.constants(2)*eye(size(Tap'*Tap,1));
    Ed                             = options.constants(3)*eye(size(Tad'*Tad,1));
    Ec                             = options.constants(1)*eye(nwx+2*nx);
    
    % Construct LMI
    prob                           = iqclmi(prob,Vx,-1,Tap'*Tcp*Tap+Ep,Top*Tap);
    prob                           = iqclmi(prob,Vy, 1,Tad'*Tcd*Tad+Ed,Tod*Tad);
    prob                           = iqclmi(prob,Vc, 1,Tcc+Ec,Toc);
    prob                           = iqclmi(prob,Vxy,1,-fOblkdiag(eye(nx)));
    prob                           = iqclmi(prob,gamma,1);

    % solve LPVsyn problem
    disp(' ');
    disp('Improve the conditioning of the coupling condition...');
    prob                           = iqcsolve(prob,gamma);
    gam_cond                       = prob.gamma;

    if gam_cond == -1
         error('Suboptimal (3nd) LMIs are not feasible.');
    else
       disp(['The subobtial synthesis problem is feasible with gamma = ',num2str(ga)]) 
    end

    % call variables
    X                              = iqcdec2mat(prob,X);
    Y                              = iqcdec2mat(prob,Y);

    disp(' ');
    disp(['The norm of X  is ', num2str(norm(X),'%10.4e\n')]);
    disp(['The norm of Y  is ', num2str(norm(Y),'%10.4e\n')]);
end
disp('--------------------------------------------------------------------')

Xs.X                               = X;
Xs.dimX                            = [nwx,nx];

Ys.Y                               = Y;
Ys.dimY                            = [nwx,nx];

% Construct the extended Lyapunov matrix
Xe                                 = fLyapext(Xs,Ys,14);

% Construct the controller
Xet                                = [Xe,blkdiag(CW',zeros(nwx+2*nx,nwi));zeros(nwi,2*(nwx+nx)),DW',Mga];
T1                                 = (Xet*[blkdiag(zeros(nwx,nwx+nx),B1);blkdiag(eye(nwx+nx),E);zeros(nwi,nwx+nx+nu)])';
T2                                 = [blkdiag(zeros(nwx+nx,nwx),C1),blkdiag(eye(nwx+nx),F)];
T3                                 = fHe(Xet*[blkdiag(AW,A,zeros(nwx+nx)),[BW;B;zeros(nwx+nx,nwi)];zeros(nwo,nwx),C,zeros(nwo,nwx+nx),D;fJt(nwi,2*(nwx+nx))']);

Kt                                 = fInvproj(T1,T2,T3);
K                                  = (eye(size(Kt,1))+Kt*blkdiag(zeros(nwx+nx),Dyu))\Kt;
K                                  = ss(K(1:nwx+nx,1:nwx+nx),K(1:nwx+nx,nwx+nx+1:end),K(nwx+nx+1:end,1:nwx+nx),K(nwx+nx+1:end,nwx+nx+1:end));
D                                  = fCleanupsys(ss(K.d),1e-14);
K.d                                = D.d;
end