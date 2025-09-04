function [L,ga] = fAWsyn(Pl,nout,nin,options)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        25-03-2021
% 
% -------------------------------------------------------------------------
%
% Description: This function performs an anti-windup compensator synthesis
%
% Syntax:      [K,gamma]    = fHinfsyn(Pl,nout,nin)
%              [K,gamma]    = fHinfsyn(Pl,nout,nin,options)
%
% Usage:       This function solves a particular version of an anti-windup
%              compensator design problem as reported in [1].
%
%              [1] E.F. Mulder, M.V. Kothare, M. Morari, "Multivariable
%              anti-windup controller synthesis using linear matrix 
%              inequalities", Automatica 37 (2001) 1407-1416.
%
%              Given the closed loop generalized plant interconnection
%
%                            |-----------|
%                  z <-------|           |<--------- w
%                            |     P     |
%                        |---|           |<--|
%                        |   |-----------|   |
%                      y |                   | u
%                        |      |-----|      |
%                        |----->|  K  |------|
%                               |-----|
%
%              where [z;y] = P [w;u] is the (weighted) generalized plant
%              and K = (Ak,Bk,Ck,Dk) the controller, this algorithm
%              designs an anti-windup compensator, L, for the extended
%              interconnection 
%
%
%                            |-----------|
%                  z <-------|           |<--------- w
%                            |     P     | 
%                        |---|           |<--O<-------|----------|
%                        |   |-----------|   ^ -      | p        |
%                        |                   |     |------|      |  y
%                        |                   |     |  dz  |      v
%                        |     |-------|     |     |------|   |-----|
%                        |---->|       |     |        ^       |  L  |
%                              |   Ke  |-----|--------| q     |-----|
%                        |---->|       |                         |
%                        |     |-------|                         |  u
%                        |                                       |
%                        |---------------------------------------|
%
%              Here Ke is extended as
%
%                     [ Ak | Bk  [I 0] ]
%                Ke = [----|-----------]
%                     [ Ck | Dk  [0 I] ]
%
%              and sat(.) is the saturation nonlinearity, which can be
%              represented as  sat(.) = 1-dz(.) where dz(.) is a deadzone
%              nonlinearity.
%
%              As input:
%
%                # The plant Pl
%
%                       {       [ A | Bp , Bw , Bu ]
%                       { [q]   [---|--------------] [p]
%                  Pl = { [z] = [Cq | Dqp, Dqw, Dqu] [w]
%                       { [y]   [Cz | Dzp, Dzw, Dzu] [u]
%                       {       [0  |   I,   0,   0]
%
%                  with 
%                    - uncertainty in/output:  p, q
%                    - performance in/output:  w, z
%                    - control in/output:      y, u
%
%                  Here the performance channel satisfies:
%
%                  int_0^\infty [*][ I,       0 ][z(t)]dt\geq0
%                               [*][ 0, -ga^2 I ][w(t)]
%
%                  for all t\geq0 and where gamma is minimized.
%
%                  On the other hand, in the case that nin and nout are
%                  defined by
%
%                  nin = [np;nw;nu]  and  nout = [nq;nz;ny]
%
%                  where:
%
%                    - np denotes the number of uncertainty inputs,
%                    - nw denotes the number of external disturbance inputs,
%                    - nu denotes the number of control inputs.
%
%                  and where
%
%                    - nq denotes the number of uncertainty outputs,
%                    - nz denotes the number of performance outputs,
%                    - ny denotes the number of measurements outputs.
%
%                # The optional structure "options", where:
%
%                   - options.alpha       If specified, this specifies the
%                                         sector constraint [0,aplha]
%                                         (default = 1)
%
%                   - options.subopt      If larger than 1, the algorithm
%                                         computes a suboptimal solution
%                                         options.subopt*gamma while
%                                         minimizing the norms on the
%                                         LMI variables (default = 1.01).
%
%                   - options.constants   options.constants = [f1,f2,f3]
%                                         is a verctor whose elements
%                                         perturb the LMIs with LMIi < -fiI
%                                         (Thus fi should be small)
%                                         (default = [0,0,0]).
%
%                   - options.Parser      Use LMIlab or Yalmip to solve the
%                                         optimization problem (default =
%                                         'LMIlab') 
%
%                   - options.Solver      Use mincx in case LMIlab is
%                                         considered or another solver
%                                         (e.g. sedumi, sdpt3, mosek, etc.)
%                                         in case Yalmip is used.
%
%                   - options.FeasbRad    Bound variables (Default = 1e6)
%
%                   - options.Terminate   can be used to change the
%                                         LMIlab-solver details (see help:
%                                         minxc) (default = 20) 
%
%                   - options.RelAcc      can be used to change de
%                                         LMIlab-solver details (see help:
%                                         minxc) (default = 1e-4)
%
%              As Output:
%
%                # The anti-windup compensator K, which stabilizes the
%                  closed-loop system lft(Pl,K) while the following
%                  performance criteria is satisfied:
%
%                  int_0^\infty [z(t)]^T[ I,       0 ][z(t)]dt\geq0
%                               [w(t)]  [ 0, -ga^2 I ][w(t)]
%
%                  for all t\geq0.
%
%                # The L2-infinity norm "ga"
%
% -------------------------------------------------------------------------
if nargin == 3
    options.alpha                  = 1;
    options.subopt                 = 1.01;    
    options.constants              = [0,0,0];    
    options.Parser                 = 'LMIlab';
    options.Solver                 = 'mincx';
    options.FeasbRad               = 1e6;
    options.Terminate              = 20;
    options.RelAcc                 = 1e-4;
elseif nargin == 4
    if isfield(options,'alpha')
        if options.alpha > 0
            options.alpha          = options.alpha;
        else
            options.alpha          = 1;
        end
    else
        options.alpha              = 1;
    end
    if isfield(options,'subopt')
        if options.subopt >= 1
            options.subopt         = options.subopt;
        else
            options.subopt         = 1.01;
        end
    else
        options.subopt             = 1.01;
    end
    if isfield(options,'constants')
        if min(options.constants) < 0
            options.constants      = [0,0,0];
        else
            options.constants      = options.constants;
        end
    else
        options.constants          = [0,0,0];
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
    if isfield(options,'FeasbRad')
        if options.FeasbRad > 0
            options.FeasbRad       = options.FeasbRad;
        else
            options.FeasbRad       = 1e6;
        end
    else
        options.FeasbRad           = 1e6;
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

% Define dimensional info
nx                                 = size(Pl.a,1);
np                                 = nin(1);
nw                                 = nin(2);
nu                                 = nin(3);
nq                                 = nout(1);
nz                                 = nout(2);
ny                                 = nout(3);

if np ~= nq
    error('Error: The IO dimensions of the uncertainty channel should be equal');
end

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

if norm(Cy) ~= 0 || norm(Dyw) ~= 0 || norm(Dyu) ~= 0
    error('The state space matrices Cy, Dyw, and Dyu should be zero');
end
if norm(Dyp) ~= 1 && trace(Dyp) ~=np
   error('The state space matrix Dyp should be the identity matrix'); 
end

% Define sector constraints [0,alp]
alp                                = options.alpha;

% Define LMI problem
probOpt.FeasbRad                   = options.FeasbRad;
probOpt.Parser                     = options.Parser;
probOpt.Solver                     = options.Solver;
probOpt.Terminate                  = options.Terminate;
probOpt.RelAcc                     = options.RelAcc;
prob                               = iqcprob(probOpt);

% Define LMI variables
X                                  = iqcvar(prob,[nx,nx],'symmetric');
Y                                  = iqcvar(prob,[nu,np],'full');
M                                  = diag(iqcvar(prob,[1,np],'full'));
gamma                              = iqcvar(prob,[1,1],'symmetric');
delta                              = iqcvar(prob,[1,1],'symmetric');

At                                 = [A,alp*Bp,Bu;zeros(nw,nx+np+nu);...
                                     Cq,alp*Dqp-eye(nq),Dqu;Cz,alp*Dzp,Dzu;...
                                     zeros(nq,nx),alp*eye(np),zeros(nq,nu)];
Xt                                 = blkdiag(X,[zeros(np+nu,nw),[M;Y],zeros(np+nu,nz+nq)]);

Mt                                 = blkdiag(-gamma*eye(nw+nz),-delta*eye(np));
Bt                                 = fHe([zeros(nx+nw+np+nz+nq,nx),...
                                         [Bw;zeros(nw);Dqw;Dzw;zeros(nq,nw)],...
                                         zeros(nx+nw+np+nz+nq,np+nz+nq)]);

Bt_eps                             = Bt + options.constants(1)*eye(size(Bt,2));

Jt1                                = eye(size(At,1));
Jt2                                = [zeros(nw,nx),eye(nw),zeros(nw,np+nz+nq);...
                                      zeros(nz+np,nx+nw+np),eye(nz+np)];
Tt                                 = [Jt1;At';Jt2];
Wt                                 = blkdiag(oblkdiag(Xt'),Mt);

% Construct LMI
prob                               = iqclmi(prob,Wt,-1,-Bt_eps,Tt);
prob                               = iqclmi(prob,X,1,options.constants(2)*eye(nx));
prob                               = iqclmi(prob,M,1,options.constants(3)*eye(np));
prob                               = iqcsolve(prob,gamma);

ga                                 = prob.gamma;

if ga == -1
    error('LMIs are not feasible.');
else
    disp(['A feasible solution was found with gamma = ',num2str(ga)])
end

if options.subopt > 1
    % Define LMI problem
    probOpt.FeasbRad               = options.FeasbRad;
    probOpt.Parser                 = options.Parser;
    probOpt.Solver                 = options.Solver;
    probOpt.Terminate              = options.Terminate;
    probOpt.RelAcc                 = options.RelAcc;
    prob                           = iqcprob(probOpt);
    
    % set suboptimal gamma
    ga                             = options.subopt*ga;
    
    % define LMI variable
    X                              = iqcvar(prob,[nx,nx],'symmetric');
    Y                              = iqcvar(prob,[nu,np],'full');
    M                              = diag(iqcvar(prob,[1,np],'full'));
    gamma                          = iqcvar(prob,[1,1],'symmetric');
    delta                          = iqcvar(prob,[1,1],'symmetric');
     
    % Bound LMIvariables
    prob                           = fBoundvar(prob,Y,gamma);
    prob                           = fBoundvar(prob,M,gamma);
    
    At                             = [A,alp*Bp,Bu;zeros(nw,nx+np+nu);...
                                      Cq,alp*Dqp-eye(nq),Dqu;Cz,alp*Dzp,Dzu;...
                                      zeros(nq,nx),alp*eye(np),zeros(nq,nu)];
    Xt                             = blkdiag(X,[zeros(np+nu,nw),[M;Y],zeros(np+nu,nz+nq)]);

    Mt1                            = -delta*eye(np);
    Mt2                            = -ga*eye(nw+nz);
    Bt1                            = fHe([zeros(nx+nw+np+nz+nq,nx),...
                                         [Bw;zeros(nw);Dqw;Dzw;zeros(nq,nw)],...
                                         zeros(nx+nw+np+nz+nq,np+nz+nq)]);
    Bt2                            = [zeros(nw,nx),eye(nw),zeros(nw,np+nz+nq);zeros(nz,nx+nw+np),eye(nz),zeros(nz,nq)];
                                  
    Bt_eps                         = Bt1 + Bt2'*Mt2*Bt2 + options.constants(1)*eye(size(Bt,2));

    Jt1                            = eye(size(At,1));
    Jt2                            = [zeros(np,nx+nw+np+nz),eye(np)];
    Tt                             = [Jt1;At';Jt2];
    Wt                             = blkdiag(oblkdiag(Xt'),Mt1);

    % Construct LMI
    prob                           = iqclmi(prob,Wt,-1,-Bt_eps,Tt);
    prob                           = iqclmi(prob,X,1,options.constants(2)*eye(nx));
    prob                           = iqclmi(prob,M,1,options.constants(3)*eye(np));
    prob                           = iqclmi(prob,gamma,1);
    prob                           = iqclmi(prob,delta,1);

    % solve Hinfsyn problem
    prob                           = iqcsolve(prob,gamma);
    ga_bnd                         = prob.gamma;
    if ga_bnd == -1
        error('Suboptimal (2nd) LMIs are not feasible!');
    else
        disp('The suboptimal syntehsis problem is feasible');
    end
end

% call variables
Y                                  = iqcdec2mat(prob,Y);
M                                  = iqcdec2mat(prob,M);
L                                  = Y/M;
end