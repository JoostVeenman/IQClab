function [Nrefined,RefinedBounds] = fRefine(N,sigs,params,options)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        15-08-2025
% 
% -------------------------------------------------------------------------
%
% Description: This function performs a data-based paramter refinedment
%              based on the algorithm described in the paper: 
%
%              Tobias Holicki, Carsten W. Scherer. Data-Based Refinement of
%              Parametric Uncertainty Descriptions. Authorea. December 16,
%              2024. DOI: 10.22541/au.173434872.23299893/v1
%
%              The aim of this function is to refine uncertain LTI
%              parameter bounds by means of using measurement data taken
%              from the real system.
%
% Syntax:      [Nrefined,RefinedBounds] = fRefine(N,sigs,params)
%              [Nrefined,RefinedBounds] = fRefine(N,sigs,params,options)
%
% Usage:       [Nrefined,RefinedBounds] = fRefine(N,sigs,params) performs a
%              parameter refinement for the uncertain system N. Here:
%
%                1.) N is an LFT [Md,Delta] = lftdata(N) where:
%
%                      - Md is the open-loop LTI discrete time plant model 
%                      - Delta is the uncertainty block consisting of: 
%                              - LFT parametric uncertainties whose bounds
%                                are to be refined,
%                              - LTV (or LTI) parametric uncertainties
%                                whose bounds are NOT to be refined.
%
%                2.) sigs is a structure containing the signals:
%
%                      - sigs.r: The input vector r = col(r(1),...,r(n)).
%                      - sigs.x: The state vector x = col(x(1),...,x(n)).
%                                (only required for Case = 'ISdata', see
%                                below for details).
%                      - sigs.y: The output vector y = col(y(1),...,y(n))
%                                (only required for Case = 'IOdata', see
%                                below for details).
%                      - sigs.n: The maximum noise magnitude |noise(t)|< n
%                                for all t > 0 (note: this is a scaler).
%
%                3.) params is a cell with the LTI parameter uncertainties
%                    that need to be refined.
%
%                      - params = {'de1',...,'deN'}
%
%                When specifying the 4th input, you can specify your
%                preferences under which the refinement shall be performed.
%
%                4.) The following options can be specified:
%
%                      - options.Case      : Specify the algorithm to be
%                                            used:
%                                             1.) 'IOdata' to refine
%                                                  parameter bounds using
%                                                  input/output data
%                                                  correpsonding to Thm 2,3
%                                             2.) 'ISdata' to refine
%                                                  parameter bounds using
%                                                  input/state data
%                                                  corresponding to Thm 1.
%
%                      - options.DGstruct  : Specify the underlying
%                                            structure of the DG
%                                            multiplier (more structure for
%                                            faster computations at the
%                                            price of potential
%                                            conservatism: 
%                                             1.) 'full' complete DG-
%                                                  scalings.
%                                             2.) 'struct' structured DG-
%                                                  scalings (see Rmk 5 of
%                                                  paper), default.
%                                             3.) 'diag' diagonal DG
%                                                  scalings.
%
%                      - options.nlift     : The lifting horizon nlift\in N
%                                            being a positive integer. The
%                                            higher nlift, the higher the
%                                            computational load.
%
%                      - options.Sx        : Initial state ellipsoid shape
%                                            (only relevant for Case
%                                            'IOdata').
%
%                      - options.cx        : Initial state ellipsoid centre
%                                            (only relevant for Case
%                                           'IOdata').
%
%                      - options.Parser    : Select the parser 'LMIlab' or
%                                            'Yalmip', default is 'LMIlab'.
%
%                      - options.Solver    : Select the solver. This will
%                                            be 'mincx' in case the parser
%                                            is set to 'LMIlab' and any
%                                            SPD solver supported by
%                                            Yalmip.
%
%                      - options.RelAcc    : Solver Relative accuracy,
%                                            default is 1e-3 (see mincx for
%                                            further details).
%
%                      - options.MaxNumIter: Solver maximum number of
%                                            iterations, default is 500
%                                            (see mincx for further
%                                            details).
%
%                      - options.FeasbRad  : Solver feasibility radius,
%                                            default 1e9 (see mincx for
%                                            further details).
%
%                      - options.Terminate : Solver terminal condition,
%                                            default is 100 (see mincx for
%                                            further details).
%
%                      - options.eps       : Perturb all LMIs by small esp
%                                            (i.e., LMI_i > eps I) default
%                                            is 0.
%
%                As output, the function provides:
%
%                1.) Nrefined:      The LFT plant model with refined bounds
%                2.) RefinedBounds: The refined parameter bounds
%
% -------------------------------------------------------------------------

% Set/Check input options
if nargin == 3
    probOpt.Parser     = 'LMIlab';             % LMI parser type
    probOpt.Solver     = 'mincx';              % Solver type
    probOpt.RelAcc     = 1e-3;                 % Solver Rel. acc. (see help mincx)
    probOpt.MaxNumIter = 500;
    probOpt.FeasbRad   = 1e9;
    probOpt.Terminate  = 100;                  % Solver terminal conditions (see help mincx)
    probOpt.eps        = 0;                    % Perturb all LMIs by small eps (e.g., LMI > eps I)
    Case               = 'IOdata';             % Input case IO
    DGstruct           = 'struct';             % DG scalings structure ('full','struct', or 'diag')
    nlift              = 1;                    % Lifting horizing j\in N.
    Sx                 = eye(size(N.a,1));     % Initial state ellipsoid shape
    cx                 = zeros(size(N.a,1),1); % Initial state ellipsoid centre
elseif nargin == 4
    if isfield(options,'Parser')
        if strcmp(options.Parser,'LMIlab') || strcmp(options.Parser,'Yalmip')
            probOpt.Parser = options.Parser;
        else
            probOpt.Parser = 'LMIlab';
            disp('Error: The field "Parser" can be either "LMIlab" or "Yalmip". Field set to default value, i.e., "LMIlab".')
        end
    else
        probOpt.Parser = 'LMIlab';
    end
    if isfield(options,'Solver')
        if strcmp(options.Parser,'LMIlab')
            probOpt.Solver = 'mincx';
        elseif strcmp(options.Parser,'Yalmip')
            probOpt.Solver = options.Solver;
        end
    else
        probOpt.Parser = 'LMIlab';
        probOpt.Solver = 'mincx';
    end
    if isfield(options,'RelAcc')
        if options.RelAcc>0
            probOpt.RelAcc = options.RelAcc;
        else
            probOpt.RelAcc = 1e-3;
        end
    else
        probOpt.RelAcc = 1e-3;
    end
    if isfield(options,'MaxNumIter')
        if options.MaxNumIter>0
            probOpt.MaxNumIter = MaxNumIter;
        else
            probOpt.MaxNumIter = 500;
        end
    else
        probOpt.MaxNumIter = 500;
    end
    if isfield(options,'FeasbRad')
        if options.FeasbRad > 0
            probOpt.FeasbRad = options.FeasbRad;
        else
            probOpt.FeasbRad = 1e9;
        end
    else
        probOpt.FeasbRad = 1e9;
    end
    if isfield(options,'Terminate')
        if options.Terminate > 0
            probOpt.Terminate = options.Terminate;
        else
            probOpt.Terminate = 100;
        end
    else
        probOpt.Terminate = 100;
    end
    if isfield(options,'eps')
        if options.eps > 0
            probOpt.eps = options.eps;
        else
            probOpt.eps = 0;
        end
    else
        probOpt.eps = 0;
    end
    if isfield(options,'Case')
        if strcmp(options.Case,'ISdata') || strcmp(options.Case,'IOdata')
            Case = options.Case;
        else
            Case = 'IOdata';
            disp('Error: The field "Case" can be either IOdata or ISdata. Field set to default value, i.e., "IOdata".')
        end
    else
        Case = 'IOdata';
    end
    if isfield(options,'LiftingNo')
        if isreal(options.LiftingNo) && options.LiftingNo >=1 && floor(options.LiftingNo) == options.LiftingNo
            nlift = options.LiftingNo;
        else
            nlift = 1;
            disp('Error: The field "LiftingNo" must be a positive integer. Field set to default value, i.e., 1.');
        end
    else
        nlift = 1;
    end
    if isfield(options,'DGstruct')
        if strcmp(options.DGstruct,'full') || strcmp(options.DGstruct,'struct') || strcmp(options.DGstruct,'diag')
            DGstruct = options.DGstruct;
        else
            DGstruct = 'struct';
            disp('Error: The field "DGstruct" can be either "full", "struct", or "diag". Field set to default value, i.e., "struct".');
        end
    else
        DGstruct = 'struct';
    end
    if isfield(options,'InitEll_Sx')
        if isreal(options.InitEll_Sx) && size(options.InitEll_Sx,1) == size(options.InitEll_Sx,2) && ...
                size(options.InitEll_Sx,1) == size(N.a,1) && ...
                norm(options.InitEll_Sx - diag(diag(options.InitEll_Sx))) == 0 ...
                && all(diag(options.InitEll_Sx) > 0)
            Sx = options.InitEll_Sx;
        else
            Sx = eye(size(N.a,1));
        end
    else
        Sx = eye(size(N.a,1));
    end
    if isfield(options,'InitEll_cx')
        if isreal(options.InitEll_cx) && size(options.InitEll_cx,2) == 1 && size(options.InitEll_Sx,1) == size(N.a,1)
            cx = options.InitEll_cx;
        else
            cx = zeros(size(N.a,1),1);
        end
    else
        cx = zeros(size(N.a,1),1);
    end
end

% Check validity of the structure sigs with the signals rt, yt, xt, nt
if ~isfield(sigs,'rt')
    error('To run the parameter refinement algorithm you need to provide the input trajectory sigs.rt');
else
    r                  = sigs.rt;
end
if ~isfield(sigs,'nt')
    error('To run the parameter refinement algorithm you need to provide the maximum noise magnitude sigs.nt');
else
    nt                 = sigs.nt;
end
if strcmp(Case,'ISdata')
    if ~isfield(sigs,'xt')
        error('To run the ISdata case you need to provide the state trakectory sigs.xt');
    else
        x              = sigs.xt;
    end
elseif strcmp(Case,'IOdata')
    if ~isfield(sigs,'yt')
        error('To run the IOdata case you need to provide the output trajectory sigs.yt');
    else
        y              = sigs.yt;
    end
end
%-------------------------------------------------------------------------%
%                         Perform preparations                            %
%-------------------------------------------------------------------------%

% Rearrange the uncertainty channels according to 'params'
[M,d]                  = fLFTdata(N);

% Check if N is a discrete time system
if M.Ts == 0
    error('N must be a discrete time system');
end

for i = 1:length(params)
    k(i,:)             = contains(d.name,params{i});
    n                  = find(k(i,:));
    e.name{i}          = d.name{n};
    e.occurrences(i)   = d.occurrences(n);
    e.bounds(i,1:2)    = d.bounds{n};
    e.in{i}            = d.in{n};
    e.out{i}           = d.out{n};
end

kn                     = abs(null(double(k)))';

for j = 1:size(kn,1)
    n                  = find(kn(j,:));
    e.name{i+j}        = d.name{n};
    e.occurrences(i+j) = d.occurrences(n);
    e.bounds(i+j,1:2)  = d.bounds{n};
    e.in{i+j}          = d.in{n};
    e.out{i+j}         = d.out{n};
end

f                      = e.in;
nf                     = length(f);
f{nf+1}                = d.pin;
g                      = e.out;
ng                     = length(g);
g{ng+1}                = d.pout;

[mo,mi]                = size(M);
To                     = [];
Ti                     = [];
for i = 1:length(f)
    [Toi,Tii]          = fT(mo,mi,{g{i}},{f{i}});
    To                 = [To;Toi];
    Ti                 = [Ti,Tii];
end

if sum(sum(To) > 1) || sum(sum(To.') > 1) || sum(sum(Ti) > 1) || sum(sum(Ti.') > 1)
   error('Error: At least 1 plant in- or output channel was assigned twice or more.');
end
Mr                     = To*M*Ti;
h                      = e;

ijk                    = 1;
dfg                    = 1;
for i = 1:length(h.name)
    h.in{i}            = ijk:ijk + h.occurrences(i)-1;
    h.out{i}           = dfg:dfg + h.occurrences(i)-1;
    ijk                = ijk + h.occurrences(i);
    dfg                = dfg + h.occurrences(i);
end
horg                   = h;                         % save original structure h

% Define dimensions
n                      = length(params);            % no. of LTI paramters
m                      = length(h.name);            % no. of LTI + LTV paramters
nw                     = sum(h.occurrences(1:n));   % no. of LTI param. unc. input channels 
np                     = sum(h.occurrences(n+1:m)); % no. of LTV param. unc. input channels
nr                     = size(sigs.rt,2);           % no. of reference input channels
nn                     = size(N.NominalValue,2)-nr; % no. of noise input channels
nz                     = sum(h.occurrences(1:n));   % no. of LTI param. unc. output channels
nq                     = sum(h.occurrences(n+1:m)); % no. of LTV param. unc. output channels
ny                     = size(N.NominalValue,1);    % no. of measurement output channels
nx                     = size(Mr.a,1);              % no. of states
niter                  = size(sigs.rt,1)-1;         % no. of iterations

% Extract system matrices
g                      = Sys2Mat(Mr,[nw,np,nn,nr],[nz,nq,ny]);

% Check if it is possible to run special IOCase when [Dyw,Dyp,Dyn]=0
flg_special_case       = 0;
if norm([g.Dyw,g.Dyp,g.Dyn]) == 0 && nlift == 1
    flg_special_case   = 1;
end

if n < m
    ijk                = 1;
    H                  = cell(1,length(h.occurrences(n+1:m)));
    for j = n+1:m
        H0             = zeros(sum(h.occurrences(n+1:m)));
        nh0            = h.in{j}-nw;
        H0(nh0,nh0)    = eye(h.occurrences(j));
        H{ijk}         = H0;
        ijk            = ijk + 1;
    end
    La                 = polydec(pvec('box',h.bounds(n+1:m,:)))';
    De                 = cell(1,size(La,1));
    for i = 1:size(La,1)
        De{i}          = zeros(size(H{1}));
        for j = 1:length(H)
            De{i}      = De{i} + La(i,j)*H{j};
        end
    end
elseif n == m  
    De                 = {};
end

% Lift the system if LiftingNo > 1
nno                    = nn;                        % store old nn
if nlift > 1
    % Lift system matrices
    g                  = LiftSysTot(g,nlift,h.occurrences(1:n));

    % Update dimensions
    nw                 = size(g.Bw,2);              % no. of LTI param. unc. input channels 
    np                 = size(g.Bp,2);              % no. of LTV param. unc. input channels
    nn                 = size(g.Bn,2);              % no. of noise input channels
    nr                 = size(g.Br,2);              % no. of reference input channels
    nz                 = size(g.Cz,1);              % no. of LTI param. unc. output channels
    nq                 = size(g.Cq,1);              % no. of LTV param. unc. output channels
    ny                 = size(g.Cy,1);              % no. of measurement output channels
    niter              = niter-nlift;               % no. of iterations

    % Update uncertainty information in the structure h
    h.occurrences      = h.occurrences*nlift;
    ncntr              = 1;
    for i = 1:length(h.name)
        h.in{i}        = ncntr:ncntr+h.occurrences(i)-1;
        h.out{i}       = ncntr:ncntr+h.occurrences(i)-1;
        ncntr          = ncntr + h.occurrences(i);
    end
end

RefinedBounds          = cell(1,niter);             % Set output cell with Refined parameter bounds

switch Case
case 'ISdata'
%-------------------------------------------------------------------------%
%                   Theorem 1 of [Holicki,2025]                           %
%-------------------------------------------------------------------------%

% Step 1: Initialization (set t <- 0, a <- a0, b <- b0 and choose an
% ellipsoid E(S,C) containing B(a,b) tightly by using Corollary 1 of
% [Holicki,2025].
[S,c]                  = CoverBoxWithEllipse(h.bounds(1:n,1),h.bounds(1:n,2));
Scal                   = SelectionMatrix(nw,n,h);

for i = 1:niter
    updateCounter(i,niter);

    % Step 2: Find a smallest ellipsoid E(Snew,cnew) that is guaranteed to
    % contain the uncertainties delta_* based on all the available information
    % at time t including \delta_*\in E(Snew,cnew)\cap B(anew,bnew).

    prob               = iqcprob(probOpt);

    rt                 = r(i:i+nlift-1,:)'; rt  = rt(:);
    xt                 = x(i,:)';           xt  = xt(:);
    xtn                = x(i+1,:)';         xtn = xtn(:);

    R                  = [1,zeros(1,n+nw+np+nn);...
                         -c,eye(n),zeros(n,nw+np+nn);...
                         g.A*xt+g.Br*rt-xtn,zeros(nx,n),g.Bw,g.Bp,g.Bn;...
                         1,zeros(1,n+nw+np+nn);...
                         g.Cz*xt+g.Dzr*rt,zeros(nz,n),g.Dzw,g.Dzp,g.Dzn;...
                         zeros(n+nw,1),eye(n+nw),zeros(n+nw,np+nn);...
                         g.Cq*xt+g.Dqr*rt,zeros(nq,n),g.Dqw,g.Dqp,g.Dqn;...
                         zeros(np,1+n+nw),eye(np),zeros(np,nn);...
                         1,zeros(1,n+nw+np+nn);...
                         zeros(nn,1+n+nw+np),eye(nn)];
    Tsc                = DGscalingMatrix(h.bounds(1:n,:),h.occurrences(1:n)+1);
    Rcal               = blkdiag(1,chol(S^-1),eye(nx),Tsc*Scal,eye(nq+np),nt*ones(nlift,1),eye(nn))*R;

    Snew               = diag(iqcvar(prob,[n,1],'full'));
    cnew               = iqcvar(prob,[n,1],'full');
    tau_eps            = iqcvar(prob,[1,1],'full');
    tau_x              = iqcvar(prob,[1,1],'full');
    tau_n              = iqcvar(prob,[nlift,1],'full');
    gamma              = iqcvar(prob,[1,1],'full');
    [M,prob]           = DGvars(prob,h.occurrences(1:n)+1,DGstruct);
    [P,prob]           = Pvars(prob,De,nlift);

    Mcal               = blkdiag(tau_eps,-tau_eps*eye(n),tau_x*eye(nx),M,P,diag(tau_n),-kron(diag(tau_n),eye(nno)));
    V                  = blkdiag([0,-cnew';-cnew,-Snew],Mcal);
    O1                 = [blkdiag(Rcal(1,:),eye(n));Rcal,zeros(size(Rcal,1),n)];
    Co1                = [Rcal(1,:)'*Rcal(1,:),[zeros(1,n);-eye(n);zeros(nw+np+nn,n)];[zeros(1,n);-eye(n);zeros(nw+np+nn,n)]',zeros(n)];
    W                  = blkdiag(Snew,-gamma);
    O2                 = ones(Snew.Dim(1,1)+1,1);

    prob               = iqclmi(prob,tau_eps,1);
    prob               = iqclmi(prob,diag(tau_n),1);
    prob               = iqclmi(prob,gamma,1);
    prob               = iqclmi(prob,V,-1,Co1,O1);
    prob               = iqclmi(prob,W,-1,0,O2);
    prob.Display       = 'offe';
    prob               = iqcsolve(prob,gamma);

    Snew               = iqcdec2mat(prob,Snew);
    cnew               = iqcdec2mat(prob,cnew);

    % Step 3: Determine a box B(anew,bnew) contraining E(Snew,cnew)\cap
    % B(anew,bnew) based on Corollary 2 of [Holicki,2025]
    [anew,bnew]        = CoveringEllipsoidWithBox(Snew,cnew);
    
    % Step 4: Update a <- anew, b<- bnew, S <- Snew, c <- cnew
    for j = 1:n
        h.bounds(j,1)  = max(h.bounds(j,1),anew(j));
        h.bounds(j,2)  = min(h.bounds(j,2),bnew(j));
    end
 
    S                  = Snew;
    c                  = cnew;
    RefinedBounds{i}   = h.bounds(1:n,:);
    
    % Step 5: If t = niter stop. Otherwise set t <- t+1 and repeat algorithm
end

case 'IOdata'
%-------------------------------------------------------------------------%
%                   Theorem 2 and 3 of [Holicki,2025]                     %
%-------------------------------------------------------------------------%

% Step 1: Initialization, Set t <- 0, a <- a0, b <- b0 and choose an
% ellipsoid E(S,c) containing B(a,b) tightly by using Corollary 1 of
% [Holicki,2025].
[S,c]                  = CoverBoxWithEllipse(h.bounds(1:n,1),h.bounds(1:n,2));
Scal                   = SelectionMatrix(nw,n,h);

for i = 1:niter
    updateCounter(i,niter);
    
    % Step 2 (Theorem 2 of [Holicki, 2025]): Find a 'smallest' ellipsoid
    % E(S_{x,new},c_{x,new}) that is guaranteed to contain x(t+1) based on
    % all the available information including x(t)\in E(Sx,cx) and
    % delta_*\in E(S,c)\cap B(a,b). 
    rt                 = r(i:i+nlift-1,:)'; rt  = rt(:);
    rtn                = r(i+1:i+nlift,:)'; rtn = rtn(:);
    yt                 = y(i:i+nlift-1,:)'; yt  = yt(:);
    ytn                = y(i+1:i+nlift,:)'; ytn = ytn(:);

    prob               = iqcprob(probOpt);
    R1                 = [g.Br*rt,g.A,g.Bw,g.Bp,g.Bn;...
                         1,zeros(1,nx+nw+np+nn);...
                         -cx,eye(nx),zeros(nx,nw+np+nn);...
                         g.Dyr*rt-yt,g.Cy,g.Dyw,g.Dyp,g.Dyn];
    R2                 = [g.Dzr*rt,g.Cz,g.Dzw,g.Dzp,g.Dzn;...
                         zeros(nw,1+nx),eye(nw),zeros(nw,np+nn);...
                         g.Dqr*rt,g.Cq,g.Dqw,g.Dqp,g.Dqn;...
                         zeros(nq,1+nx+nw),eye(np),zeros(nq,nn);...
                         1,zeros(1,nx+nw+np+nn);...
                         zeros(nn,1+nx+nw+np),eye(nn)];
    R12                = [g.Cy*g.Br*rt+g.Dyr*rtn-ytn,g.Cy*[g.A,g.Bw,g.Bp,g.Bn]];

    if flg_special_case == 0
        R              = [R1;R2];
        scny           = 1;
    elseif flg_special_case == 1 % only for nlift == 1
        R              = [R1;R12;R2];
        scny           = 2;
    end
    
    Tsc                = DGscalingMatrix(h.bounds(1:n,:),h.occurrences(1:n));
    Rcal               = blkdiag(eye(nx),1,chol(Sx^-1),eye(scny*ny),Tsc,eye(nq+np),nt*ones(nlift,1),eye(nn))*R;
    Rcal1              = Rcal(1:nx,:);
    Rcal2              = Rcal(nx+1:end,:);

    Sxnew              = diag(iqcvar(prob,[nx,1],'full'));
    cxnew              = iqcvar(prob,[nx,1],'full');
    tau_epsx           = iqcvar(prob,[1,1],'full');
    tau_y              = iqcvar(prob,[1,1],'full');
    tau_n              = iqcvar(prob,[nlift,1],'full');
    gamma              = iqcvar(prob,[1,1],'full');
    
    [M,prob]           = DGvars(prob,h.occurrences(1:n),DGstruct);
    [P,prob]           = Pvars(prob,De,nlift);

    Mcal               = blkdiag(tau_epsx,-tau_epsx*eye(nx),tau_y*eye(scny*ny),M,P,diag(tau_n),-kron(diag(tau_n),eye(nno)));
    V                  = blkdiag([0,-cxnew';-cxnew,-Sxnew],Mcal);
    O1                 = [blkdiag([1,zeros(1,nx+nw+np+nn)],eye(nx));Rcal2,zeros(size(Rcal2,1),nx)];
    Co1                = [[1,zeros(1,nx+nw+np+nn)]'*[1,zeros(1,nx+nw+np+nn)],-Rcal1';-Rcal1,zeros(nx)];
    W                  = blkdiag(Sxnew,-gamma);
    O2                 = ones(Sxnew.Dim(1,1)+1,1);

    prob               = iqclmi(prob,tau_epsx,1);
    prob               = iqclmi(prob,diag(tau_n),1);
    prob               = iqclmi(prob,gamma,1);
    prob               = iqclmi(prob,V,-1,Co1,O1);
    prob               = iqclmi(prob,W,-1,0,O2);
    prob.Display       = 'offe';
    prob               = iqcsolve(prob,gamma);

    Sxnew              = iqcdec2mat(prob,Sxnew);
    cxnew              = iqcdec2mat(prob,cxnew);
    
    % Step 3(Theorem 2 of [Holicki, 2025]): Find a smallest ellipsoid
    % E(Snew,cnew) that is guaranteed to contain delta_* based on all the
    % available information including x(t)\in E(Sx,cx) and delta_*\in
    % E(S,c)\cap B(a,b) 

    prob               = iqcprob(probOpt);

    R1                 = [zeros(n,1),eye(n),zeros(n,nx+nw+np+nn);...
                         1,zeros(1,n+nx+nw+np+nn);...
                         -c,eye(n),zeros(n,nx+nw+np+nn);...
                         1,zeros(1,n+nx+nw+np+nn);...
                         g.Br*rt-cxnew,zeros(nx,n),g.A,g.Bw,g.Bp,g.Bn;...
                         1,zeros(1,n+nx+nw+np+nn);...
                         -cx,zeros(nx,n),eye(nx),zeros(nx,nw+np+nn);...
                         g.Dyr*rt-yt,zeros(ny,n),g.Cy,g.Dyw,g.Dyp,g.Dyn];
    R2                 = [1,zeros(1,n+nx+nw+np+nn);...
                         g.Dzr*rt,zeros(nz,n),g.Cz,g.Dzw,g.Dzp,g.Dzn;...
                         zeros(n,1),eye(n),zeros(n,nx+nw+np+nn);...
                         zeros(nw,1+n+nx),eye(nw),zeros(nw,np+nn);...
                         g.Dqr*rt,zeros(nq,n),g.Cq,g.Dqw,g.Dqp,g.Dqn;...
                         zeros(nq,1+n+nx+nw),eye(np),zeros(nq,nn);...
                         1,zeros(1,n+nx+nw+np+nn);...
                         zeros(nn,1+n+nx+nw+np),eye(nn)];
    R12                = [g.Cy*g.Br*rt+g.Dyr*rtn-ytn,zeros(ny,n),g.Cy*[g.A,g.Bw,g.Bp,g.Bn]];

    if flg_special_case == 0
        R              = [R1;R2];
        scny           = 1;
    elseif flg_special_case == 1
        R              = [R1;R12;R2];
        scny           = 2;
    end

    Tsc                = DGscalingMatrix(h.bounds(1:n,:),h.occurrences(1:n)+1);
    
    Rcal               = blkdiag(eye(n),1,chol(S^-1),1,chol(Sxnew^-1),1,chol(Sx^-1),eye(scny*ny),Tsc*Scal,eye(nq+np),nt*ones(nlift,1),eye(nn))*R;
    Rcal1              = Rcal(1:n,:);
    Rcal2              = Rcal(n+1:end,:);

    Snew               = diag(iqcvar(prob,[n,1],'full'));
    cnew               = iqcvar(prob,[n,1],'full');
    tau_del            = iqcvar(prob,[1,1],'full');
    tau_epsxn          = iqcvar(prob,[1,1],'full');
    tau_epsx           = iqcvar(prob,[1,1],'full');
    tau_y              = iqcvar(prob,[1,1],'full');
    tau_n              = iqcvar(prob,[nlift,1],'full');
    gamma              = iqcvar(prob,[1,1],'full');
    
    M                  = DGvars(prob,h.occurrences(1:n)+1,DGstruct);
    [P,prob]           = Pvars(prob,De,nlift);

    Mcal               = blkdiag(tau_del,-tau_del*eye(n),tau_epsxn,-tau_epsxn*eye(nx),tau_epsx,-tau_epsx*eye(nx),tau_y*eye(scny*ny),M,P,diag(tau_n),-kron(diag(tau_n),eye(nno)));
    V                  = blkdiag([0,-cnew';-cnew,-Snew],Mcal);
    O1                 = [blkdiag([1,zeros(1,n+nx+nw+np+nn)],eye(n));Rcal2,zeros(size(Rcal2,1),n)];
    Co1                = [[1,zeros(1,n+nx+nw+np+nn)]'*[1,zeros(1,n+nx+nw+np+nn)],-Rcal1';-Rcal1,zeros(n)];
    W                  = blkdiag(Snew,-gamma);
    O2                 = ones(Snew.Dim(1,1)+1,1);

    prob               = iqclmi(prob,tau_del,1);
    prob               = iqclmi(prob,tau_epsxn,1);
    prob               = iqclmi(prob,tau_epsx,1);
    prob               = iqclmi(prob,diag(tau_n),1);
    prob               = iqclmi(prob,gamma,1);
    prob               = iqclmi(prob,V,-1,Co1,O1);
    prob               = iqclmi(prob,W,-1,0,O2);
    prob.Display       = 'offe';
    prob               = iqcsolve(prob,gamma);

    Snew               = iqcdec2mat(prob,Snew);
    cnew               = iqcdec2mat(prob,cnew);
 
    % Step 4: Determine a box B(anew,bnew) containing E(Snew,cnew)\cap
    % B(a,b) tightly by using Corrollary 2 to conclude that \delta_*\in
    % E(Snew,cnew)\cap B(anew,bnew)\subset E(S,c)\cap B(a,b)
    [anew,bnew]        = CoveringEllipsoidWithBox(Snew,cnew);
    
    % Step 5: Update a <- anew, b<- bnew, S <- Snew, c <- cnew
    for j = 1:n
        lb             = max(h.bounds(j,1),anew(j));
        ub             = min(h.bounds(j,2),bnew(j));
        if lb < ub
            h.bounds(j,1) = lb;
            h.bounds(j,2) = ub;
        end
    end

    S                  = Snew;
    c                  = cnew;
    Sx                 = Sxnew;
    cx                 = cxnew;
    RefinedBounds{i}   = h.bounds(1:n,:);

    % Step 6: If t = niter stop. Otherwise set t <- t+1 and repeat algorithm
end
end

if nlift > 1
    nocc               = horg.occurrences;
else
    nocc               = h.occurrences;
end
Delr                   = [];
for j = 1:length(h.name)    
    delj               = ureal(h.name{j},0.5*(h.bounds(j,1)+h.bounds(j,2)),'Range',h.bounds(j,:)');
    Delr               = blkdiag(Delr,delj*eye(nocc(j)));
end
Nrefined = lft(Delr,Mr);
%-------------------------------------------------------------------------%
%                          auxiliary functions                            %
%-------------------------------------------------------------------------%
function [Sn,cn]       = CoverBoxWithEllipse(a,b)
    % Function to (tightly) cover (overbound) a box with an Ellipsoid
    d1                 = 0.5*(b-a);
    Sn                 = sum(d1)*diag(d1);
    cn                 = 0.5*(a+b);
end

function [an,bn]       = CoveringEllipsoidWithBox(S,c)
    % Function to (tightly) cover (overbound) Ellipsoid with a box
    Ssqrt              = sqrt(diag(S));
    an                 = c - Ssqrt(:);
    bn                 = c + Ssqrt(:);
end

function Tsc           = DGscalingMatrix(bnds,nocc)
    % Function to compute the DG scaling outer factor matrix
    Ta                 = coder.nullcopy(zeros(sum(nocc)));
    Tb                 = coder.nullcopy(zeros(sum(nocc)));
    b                  = cumsum(nocc);
    a                  = b - nocc + 1;
    for q = 1:length(nocc)
        lmn            = a(q):b(q);
        Ta(lmn,lmn)    = bnds(q,2)*eye(nocc(q));
        Tb(lmn,lmn)    = -bnds(q,1)*eye(nocc(q));
    end
    Tsc                = [Ta,-eye(sum(nocc));Tb,eye(sum(nocc))];
end

function [M,prob]      = DGvars(prob,nocc,opt)
    % Function to compute the DG IQC variables and define the constraint 
    % Mi + Mi > 0
    M                  = iqcvar(prob,[0,0],'full');
    for q = 1:length(nocc)
        if strcmp(opt,'full') % Mi is full
            Mi         = iqcvar(prob,[nocc(q),nocc(q)],'full');
        elseif strcmp(opt,'struct') % Mi is strutured (see Remark 5 of [Holicki, 2025]
            Mi         = [iqcvar(prob,[1,1],'full'),iqcvar(prob,[1,nocc(q)-1],'full');...
                         iqcvar(prob,[nocc(q)-1,1],'full'),diag(iqcvar(prob,[nocc(q)-1,1],'full'))];
        elseif strcmp(opt,'diag') % Mi is diagonal;
            Mi         = diag(iqcvar(prob,[nocc(q),1],'full'));
        end
        prob           = iqclmi(prob,oblkdiag(Mi'),1,0,[eye(nocc(q));eye(nocc(q))]);
        M              = blkdiag(M,Mi);
    end
    M                  = oblkdiag(M');
end

function [P,prob]      = Pvars(prob,De,nlift)
    % Function to compute the full block multiplier P and constraints P22<0
    % and (x)^TP([Delta^j])>0
    if isempty(De)
        P              = iqcvar(prob,[0,0],'symmetric');
    else
        [nin,nout]     = size(De{1});
        P              = iqcvar(prob,[nin+nout,nin+nout],'symmetric');
        for q = 1:length(De)
            prob       = iqclmi(prob,P,1,0,[eye(nin);De{q}]);
        end
            prob       = iqclmi(prob,P,-1,0,[zeros(nin,nout);eye(nout)]);
        if nlift > 1
            P11        = P(1:nin,1:nin);
            P12        = P(1:nin,nin+1:nin+nout);
            P22        = P(nin+1:nin+nout,nin+1:nin+nout);
            for q = 1:nlift-1
                Pi     = iqcvar(prob,[nin+nout,nin+nout],'symmetric');
                for p = 1:length(De)
                  prob = iqclmi(prob,Pi,1,0,[eye(nin);De{p}]);
                end
                prob   = iqclmi(prob,Pi,-1,0,[zeros(nin,nout);eye(nout)]);

                P11    = blkdiag(P11,Pi(1:nin,1:nin));
                P12    = blkdiag(P12,Pi(1:nin,nin+1:nin+nout));
                P22    = blkdiag(P22,Pi(nin+1:nin+nout,nin+1:nin+nout));
            end
            P          = [P11,P12;P12',P22];
        end
    end    
end

function Scal          = SelectionMatrix(nw,n,h)
    % Function to compute the selection matrix Scal
    J1                 = [];
    J2                 = [];

    E                  = eye(n);
    Inz                = eye(nw);
    
    for q = 1:n
        J1             = [J1;blkdiag(1,Inz(h.in{q},:))];
        J2             = [J2;blkdiag(E(:,q)',Inz(h.in{q},:))];
    end
    Scal               = blkdiag(J1,J2);
end

function g             = Sys2Mat(M,din,dout)
    mw                 = din(1);
    mp                 = din(2);
    mn                 = din(3);
    mr                 = din(4);
    mz                 = dout(1);
    mq                 = dout(2);
    my                 = dout(3);

    g.A                = M.a;
    g.Bw               = M.b(:,1:mw);
    g.Bp               = M.b(:,mw+1:mw+mp);
    g.Bn               = M.b(:,mw+mp+1:mw+mp+mn);
    g.Br               = M.b(:,mw+mp+mn+1:mw+mp+mn+mr);

    g.Cz               = M.c(1:mz,:);
    g.Cq               = M.c(mz+1:mz+mq,:); 
    g.Cy               = M.c(mz+mq+1:mz+mq+my,:);

    g.Dzw              = M.d(1:mz,1:mw);
    g.Dzp              = M.d(1:mz,mw+1:mw+mp);
    g.Dzn              = M.d(1:mz,mw+mp+1:mw+mp+mn);
    g.Dzr              = M.d(1:mz,mw+mp+mn+1:mw+mp+mn+mr);

    g.Dqw              = M.d(mz+1:mz+mq,1:mw);
    g.Dqp              = M.d(mz+1:mz+mq,mw+1:mw+mp);
    g.Dqn              = M.d(mz+1:mz+mq,mw+mp+1:mw+mp+mn);
    g.Dqr              = M.d(mz+1:mz+mq,mw+mp+mn+1:mw+mp+mn+mr);

    g.Dyw              = M.d(mz+mq+1:mz+mq+my,1:mw);
    g.Dyp              = M.d(mz+mq+1:mz+mq+my,mw+1:mw+mp);
    g.Dyn              = M.d(mz+mq+1:mz+mq+my,mw+mp+1:mw+mp+mn);
    g.Dyr              = M.d(mz+mq+1:mz+mq+my,mw+mp+mn+1:mw+mp+mn+mr);
end

function Bl            = LiftB(B,n)
    Bl                 = kron([1,zeros(1,n-1)],B);
end
function Cl            = LiftC(C,A,n)
    Cl                 = kron(ones(n,1),C);
    if n > 1
        nc             = size(C,1);
        for q = 2:n
            Cl((q-1)*nc+1:q*nc,:) = C*A^(q-1);
        end
    end
end
function Dl            = LiftD(D,C,B,A,n)
    Dl                 = kron(eye(n),D);
    if n > 1
        for q = 2:n
            Dl         = Dl + kron(diag(ones(1,n-q+1),-q+1),C*A^(q-2)*B);
        end
    end
end
function T             = Tperm(nocc,nl)
    ne                 = cumsum(nocc);
    nb                 = [1,ne(1:end-1)+1];
    T                  = zeros(sum(nocc)*nl);
    nc                 = 1;
    for q = 1:length(nocc)
        a              = ne(q)-nb(q)+1;
        b              = sum(ne(q+1:end)-nb(q+1:end)+1);
        c              = sum(ne(1:q-1)-nb(1:q-1)+1);
        Jdel           = fJ(a,b,c)';
        J              = kron(eye(nl),Jdel);
        nEc            = size(J,1);
        T(nc:nc+nEc-1,:) = J;
        nc             = nc + nEc;
    end
    T                  = T';
end

function gl            = LiftSysTot(g,n,nocc)
    % This function lifts the system matrices (A,B,C,D) as
    %
    % [*  Bi] = [*          Bi           0    ...   0   ]
    % [Ci Di]   [Ci         Di           0    ...   0   ]
    %           [CiAi       CiBi         Di   ...   0   ]
    %           [...        ...          ...  ...   ... ]
    %           [CiAi^(n-1) CiAi^(n-2)Bi ...  CiBi  Di  ]
    %
    % and permutes the LTI parametric uncertainty block (i.e.,
    % kron(Ih,Delta)) as well as the corresponding plant uncertainty
    % channels such that Delta_h = diag(delta_1 Ihn1, ... , delta_v Ihnv)
    Tp                 = Tperm(nocc,n);

    gl.A               = g.A;
    
    gl.Bw              = LiftB(g.Bw,n)*Tp;
    gl.Bp              = LiftB(g.Bp,n);
    gl.Bn              = LiftB(g.Bn,n);
    gl.Br              = LiftB(g.Br,n);
    
    gl.Cz              = Tp'*LiftC(g.Cz,g.A,n);
    gl.Cq              = LiftC(g.Cq,g.A,n);
    gl.Cy              = LiftC(g.Cy,g.A,n);
    
    gl.Dzw             = Tp'*LiftD(g.Dzw,g.Cz,g.Bw,g.A,n)*Tp;
    gl.Dzp             = Tp'*LiftD(g.Dzp,g.Cz,g.Bp,g.A,n);
    gl.Dzn             = Tp'*LiftD(g.Dzn,g.Cz,g.Bn,g.A,n);
    gl.Dzr             = Tp'*LiftD(g.Dzr,g.Cz,g.Br,g.A,n);
    
    gl.Dqw             = LiftD(g.Dqw,g.Cq,g.Bw,g.A,n)*Tp;
    gl.Dqp             = LiftD(g.Dqp,g.Cq,g.Bp,g.A,n);
    gl.Dqn             = LiftD(g.Dqn,g.Cq,g.Bn,g.A,n);
    gl.Dqr             = LiftD(g.Dqr,g.Cq,g.Br,g.A,n);
    
    gl.Dyw             = LiftD(g.Dyw,g.Cy,g.Bw,g.A,n)*Tp;
    gl.Dyp             = LiftD(g.Dyp,g.Cy,g.Bp,g.A,n);
    gl.Dyn             = LiftD(g.Dyn,g.Cy,g.Bn,g.A,n);
    gl.Dyr             = LiftD(g.Dyr,g.Cy,g.Br,g.A,n);
end

function                 updateCounter(i, niter)
    %updateCounter displays "i of niter" on a single refreshing line
    
    persistent prevStr
    
    % Build the new string
    str = sprintf('Iteration: %d of %d', i, niter);
    
    % If not the first call, erase the previous string
    if ~isempty(prevStr)
        fprintf(repmat('\b', 1, length(prevStr))); % move back
        fprintf(repmat(' ', 1, length(prevStr)));  % clear
        fprintf(repmat('\b', 1, length(prevStr))); % move back again
    end
    
    % Print the new string
    fprintf('%s', str);
    
    % Store current string for next call
    prevStr = str;
    
    % At the very last iteration, reset and move to a new line
    if i == niter
        fprintf('\n');
        prevStr = []; % reset for next time function is used
    end
end

end