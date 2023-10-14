function [E,ga] = fRobest(Pl,Delta,nout,nin,options)
% -------------------------------------------------------------------------
%
% IQClab:      Version 3.03
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available for non-commercial usage under a
%              Creative Commons (Attribution-NonCommercial-NoDerivatives
%              4.0 International (CC BY-NC-ND 4.0))license: 
%              https://creativecommons.org/licenses/by-nc-nd/4.0/
%              Commercial usage is only permitted with a commercial
%              license. For further information please visit iqclab.eu
% Author:      J.Veenman
% Date:        21-08-2023
% 
% -------------------------------------------------------------------------
%
% Description: This function performs an IQC based robust estimator
%              synthesis.
%
% Syntax:      [E,gamma]         = fRobest(Pl,Delta,nout,nin)
%              [E,gamma]         = fRobest(Pl,Delta,nout,nin,options)
%
% Usage:       This function solves the robust estimation problem
%              based on IQCs by means of an LMI solution (see iqclab.eu for
%              references.
%
%              A successful synthesis will yield an estimator, which
%              guarantees robust stability and induced L2-gain or H2
%              performance for all modelled uncertainties in the provided
%              uncertainty set.
%
%              As input one should provide:
%
%                # The (weighted) uncertain open-loop generalized plant Pl
%                  with realization
%
%                             [ A  | Bp,  Bw ]
%                       [q]   [----|---------] [p]
%                  Pl = [z] = [ Cq | Dqp, Dqw] [w]
%                       [y]   [ Cz | Dzp, Dzw]
%                             [ Cy | Dyp, Dyw]
%
%                  and with:
%                    - uncertainty input channel:  p
%                    - uncertainty output channel: q
%                    - disturbance input:          w 
%                    - performance output:         z 
%                    - measurement output:         y
%                    - A beign Hurwitz
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
%                # The input and output channel data:
%
%                  nin = [p;w]  and  nout = [q;z;y]
%
%                  where:
%
%                    - np denotes the number of uncertainty inputs
%                    - nw denotes the number of disturbance inputs
%
%                  and where
%
%                    - nq denotes the number of uncertainty outputs
%                    - nz denotes the number of performance outputs
%                    - ny denotes the number of measurements outputs
%
%                # The optional structure "options", where:
%
%                   - options.perf        Performance metric of interest:
%                                          - Induced L2-gain (default = 'L2')
%                                          - H2-norm ('H2')
% 
%                                         Note: If the option 'H2' is
%                                         selected, then Dqp must be zero
%
%                   - options.subopt      If larger than 1, the algorithm
%                                         computes a suboptimal solution
%                                         options.subopt*gamma while
%                                         minimizing the norms on the
%                                         LMI variables (Default = 1.01).
%
%                   - options.constants   options.constants = [c1,c2,c3,c4] 
%                                         is a verctor whose elements
%                                         perturb the LMIs with LMIi < -ciI
%                                         (Thus ci >= 0 should be small)
%                                         (default = 1e-9*[1,1,1,1])
%
%                                           - c1: main synthesis LMI < - c1 I 
%                                           - c2: X-Yhat > c2 I
%                                           - c3: T^T*Yhat*T - [Y,0;0,0] > c3 I
%                                           - c4: Pi11 > c4 I
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
%                # The estimator E = ss(AE,BE,CE,DE)
%
%                  For all Delta in the Delta set, this estimator 
%                  stabilizes the system while the selected performance
%                  metric is rendered less than ga.
%
%                # The worst case induced L2-gain or H_2 norm "ga".
% -------------------------------------------------------------------------
if nargin == 4
    options.perf       = 'L2';
    options.subopt     = 1.01;
    options.constants  = 1e-9*ones(1,4);
    options.Parser     = 'LMIlab';
    options.Solver     = 'mincx';
    options.FeasbRad   = 1e9;
    options.Terminate  = 0;
    options.RelAcc     = 1e-4;
    options.StrProp    = 'no';
elseif nargin == 5
    if isfield(options,'perf')
        if strcmp(options.perf,'L2') || strcmp(options.perf,'H2')
            options.perf = options.perf;
        else
            options.perf = 'L2';
        end
    else
        options.Parser = 'L2';
    end
    if isfield(options,'subopt')
        if options.subopt >= 1
            options.subopt = options.subopt;
        else
            options.subopt = 1.01;
        end
    else
        options.subopt = 1.01;
    end
    if isfield(options,'constants')
        if min(options.constants) < 0
            options.constants = 1e-9*ones(1,4);
        else
            options.constants = options.constants;
        end
    else
        options.constants = 1e-9*ones(1,4);
    end
    if isfield(options,'Parser')
        if strcmp(options.Parser,'LMIlab') || strcmp(options.Parser,'Yalmip')
            options.Parser = options.Parser;
        else
            options.Parser = 'LMIlab';
        end
    else
        options.Parser = 'LMIlab';
    end
    if isfield(options,'Solver')
        if strcmp(options.Parser,'LMIlab')
            options.Solver = 'mincx';
        elseif strcmp(options.Parser,'Yalmip')
            options.Solver = options.Solver;
        end
    else
        options.Parser = 'LMIlab';
        options.Solver = 'mincx';
    end
    if isfield(options,'FeasbRad')
        if options.FeasbRad > 0
            options.FeasbRad = options.FeasbRad;
        else
            options.FeasbRad = 1e9;
        end
    else
        options.FeasbRad = 1e9;
    end
    if isfield(options,'Terminate')
        if options.Terminate > 0
            options.Terminate = options.Terminate;
        else
            options.Terminate = 20;
        end
    else
        options.Terminate = 20;
    end
    if isfield(options,'RelAcc')
        if options.RelAcc>0
            options.RelAcc = options.RelAcc;
        else
            options.RelAcc = 1e-4;
        end
    else
        options.RelAcc = 1e-4;
    end
    if isfield(options,'StrProp')
        if strcmp(options.StrProp,'no') || strcmp(options.StrProp,'yes')
            options.StrProp = options.StrProp;
        else
            options.StrProp = 'no';
        end
    else
        options.Parser = 'no';
    end
end

% Check input vector
if abs(size(Pl.d,2) - sum(nin)) > 0 || abs(size(Pl.d,1) - sum(nout)) > 0
    error('The the vectors "nin" or "nout" are not well-defined');
end

% Define dimensional info
np                     = nin(1);
nw                     = nin(2);
nq                     = nout(1);
nz                     = nout(2);
ny                     = nout(3);

% Define LMI problem
probOpt.FeasbRad       = options.FeasbRad;
probOpt.Parser         = options.Parser;
probOpt.Solver         = options.Solver;
probOpt.Terminate      = options.Terminate;
probOpt.RelAcc         = options.RelAcc;
probOpt.eps            = 1e-9;
prob                   = iqcprob(probOpt);

if isobject(Delta)
   U{1}                = Delta;
   Delta               = U;
end

% Get sample time plant M
Ts                     = Pl.Ts;

% Permute plant in- and output uncertainty channels according to the delta
% structure 
To                     = [];
Ti                     = [];
for i = 1:length(Delta)
    [Toi,Tii]          = fT(nq,np,Delta{i}.OutputChannel,Delta{i}.InputChannel);
    To                 = [To;Toi];
    Ti                 = [Ti,Tii];
end
if sum(sum(To) > 1) || sum(sum(To.') > 1) || sum(sum(Ti) > 1) || sum(sum(Ti.') > 1)
   error('Error: At least 1 plant in- or output channel was assigned twice or more.');
end

Pl                     = blkdiag(To,eye(nz+ny))*Pl*blkdiag(Ti,eye(nw));

% Test if M is nominally stable
ePl               = eig(Pl);
if (Ts == 0) && (sum(real(ePl) >= 0)) 
    error('Error: The nominal plant is not stable.');
end

disp('-------------------------------------------------------------------');
disp(['Uncertain plant with ' num2str(size(Pl,2)) ' inputs, ' num2str(size(Pl,1)) ' outputs, and ', num2str(size(Pl.a,1)),' states.']);

% Construct IQC-multipliers
for i = 1:length(Delta)
    switch class(Delta{i})
        case 'ultis'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                prob = iqcltis_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end            
        case 'ultid'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                    set(Delta{i},'SampleTime',Ts);
                    prob = iqcltid_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'ultv'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                prob = iqcltv_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'ultv_rb'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                prob = iqcltv_rb_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'udel'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                Mu   = To'*M*Ti';
                Mu   = Mu(Delta{i}.OutputChannel{1},Delta{i}.InputChannel{1}).d;
                if norm(Mu) < 1e-14
                    set(Delta{i},'MstrictlyProp','yes');
                else
                    set(Delta{i},'MstrictlyProp','no');
                end
                set(Delta{i},'SampleTime',Ts);
                prob = iqcdel_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'usbsr'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                prob = iqcsbsr_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'usg'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                prob = iqcsg_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'ups'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                prob = iqcps_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
    end
end

% Construct extended plant compatible with the structure of "Pi"
Psi                    = get(prob,'Psi');
sc                     = get(prob,'sc');
IO                     = get(prob,'IO');

Psi1                   = sc*get(prob,'Psi1');
Psi1m                  = fss2m(Psi1);
P                      = get(prob,'P');
Y                      = iqcvar(prob,[size(Psi1.a,1),size(Psi1.a,1)],'symmetric');
kyp                    = [0,1;1,0];
V1                     = blkdiag(kron(kyp,Y),P);
prob                   = iqclmi(prob,V1,1,options.constants(4)*eye(size(Psi1m,2)),Psi1m);

% Construct extended plant compatible with the structure of "Pi"
Plex                   = fplantex(Pl,IO);

Plex                   = ss(Plex.a,Plex.b,[Plex.c(1:end-ny,:);zeros(nw,size(Plex.a,1));Plex.c(end-ny+1:end,:)],...
                            [Plex.d(1:end-ny,:);zeros(nw,size(Plex.d,2)-nw),eye(nw);Plex.d(end-ny+1:end,:)]);

Psi_sc                 = ss(sc*Psi);Psi_sc.Ts = Ts;
Psiaug                 = fAugss(Psi_sc,ss([],[],[],eye(nz+ny+nw),Ts),1);
G                      = fmultss(Psiaug,Plex);
[nGo,nGi]              = size(G);
nx                     = size(G.a,1);

A                      = G.a;
B                      = G.b;
Bp                     = G.b(:,1:end-nw);
Bw                     = G.b(:,end-nw+1:end);
npn                    = size(Bp,2);

CPsi                   = G.c(1:nGo-nz-nw-ny,:);

DPsi                   = G.d(1:nGo-nz-nw-ny,:);
DPsip                  = G.d(1:nGo-nz-nw-ny,1:end-nw);
DPsiw                  = G.d(1:nGo-nz-nw-ny,end-nw+1:end);

Cz                     = G.c(nGo-nz-nw-ny+1:1:nGo-nw-ny,:);
Dz                     = G.d(nGo-nz-nw-ny+1:1:nGo-nw-ny,:);
Dzp                    = G.d(nGo-nz-nw-ny+1:1:nGo-nw-ny,1:end-nw);
Dzw                    = G.d(nGo-nz-nw-ny+1:1:nGo-nw-ny,end-nw+1:end);

Dw                     = G.d(nGo-nw-ny+1:1:nGo-ny,:);
Cy                     = G.c(nGo-ny+1:1:end,:);
Dy                     = G.d(nGo-ny+1:1:end,:);
Dyp                    = G.d(nGo-ny+1:1:end,1:end-nw);
Dyw                    = G.d(nGo-ny+1:1:end,end-nw+1:end);

% Define LMI variables
Yhat                   = iqcvar(prob,[nx,nx],'symmetric');
X                      = iqcvar(prob,[nx,nx],'symmetric');
Ae                     = iqcvar(prob,[nx,nx],'full');
Be                     = iqcvar(prob,[nx,ny],'full');
Ce                     = iqcvar(prob,[nz,nx],'full');
if strcmp(options.StrProp,'no')
    De                 = iqcvar(prob,[nz,ny],'full');
elseif strcmp(options.StrProp,'yes')
    De                 = zeros(nz,ny);
end

gamma                  = iqcvar(prob,[1,1],'symmetric');

% Solve LMI problem
if strcmp(options.perf,'L2')
    Xtot               = [blkdiag(Yhat,[X,Ae,Be,zeros(nx,ny)]);zeros(nGi,3*nx+2*ny);zeros(nz,2*nx),Ce,De,De];
    Vtot               = blkdiag(oblkdiag(Xtot),P,-gamma*eye(nw+nz));
    OuterMain          = [eye(2*nx+nGi+nz);A,A,B,zeros(nx,nz);zeros(nx),A,B,zeros(nx,nz);...
                          eye(nx),zeros(nx,nx+nGi+nz);zeros(ny,nx),Cy,Dy,zeros(ny,nz);...
                          Cy,zeros(ny,nx+nGi+nz);CPsi,CPsi,DPsi,zeros(size(CPsi,1),nz);...
                          zeros(nw,2*nx),Dw,zeros(nw,nz);zeros(nz,2*nx+nGi),eye(nz)];
    ConstMain          = fOblkdiag([-Cz,-Cz,-Dz]')-options.constants(1)*eye(2*nx+nGi+nz);
    prob               = iqclmi(prob,Vtot,-1,ConstMain,OuterMain);
    
    Vcoupling1         = blkdiag(X,-Yhat);
    OuterCoupling1     = [eye(nx);eye(nx)];
    ConstCoupling1     = options.constants(2)*eye(nx);
    prob               = iqclmi(prob,Vcoupling1,1,ConstCoupling1,OuterCoupling1);
    
    nxpsi1             = size(prob.Psi1.a,1);
    nxpsi2             = size(prob.Psi2.a,1);
    nxPl               = size(Pl.a,1);
    Yhat11             = Yhat(1:nxpsi1,1:nxpsi1);
    Yhat13             = Yhat(1:nxpsi1,nxpsi1+nxpsi2+1:nx);
    Yhat33             = Yhat(nxpsi1+nxpsi2+1:nx,nxpsi1+nxpsi2+1:nx);
    Vcoupling2         = blkdiag([Yhat11,Yhat13;Yhat13',Yhat33],-Y);
    OuterCoupling2     = [eye(nxpsi1+nxPl);eye(nxpsi1),zeros(nxpsi1,nxPl)];
    ConstCoupling2     = options.constants(3)*eye(nxpsi1+nxPl);
    prob               = iqclmi(prob,Vcoupling2,1,ConstCoupling2,OuterCoupling2);

    disp('-------------------------------------------------------------------')
    disp('Solve main synthesis LMIs by minimizing gamma');
    prob               = iqcsolve(prob,gamma);
    ga                 = prob.gamma;

elseif strcmp(options.perf,'H2')
    if norm(DPsiw) > 1e-14
        error('Dqw must be zero'); % check if DPsiw = 0
    end
    Q                  = iqcvar(prob,[nw,nw],'symmetric');

    Xtot               = [blkdiag(Yhat,[X,Ae,Be,zeros(nx,ny)]);zeros(npn,3*nx+2*ny);zeros(nz,2*nx),Ce,De,De];
    Vtot               = blkdiag(oblkdiag(Xtot),P);
    OuterMain          = [eye(2*nx+npn+nz);A,A,Bp,zeros(nx,nz);zeros(nx),A,Bp,zeros(nx,nz);...
                          eye(nx),zeros(nx,nx+npn+nz);zeros(ny,nx),Cy,Dyp,zeros(ny,nz);...
                          Cy,zeros(ny,nx+npn+nz);CPsi,CPsi,DPsip,zeros(size(CPsi,1),nz)];
    ConstMain          = blkdiag(zeros(2*nx+npn),eye(nz))+fOblkdiag([-Cz,-Cz,-Dzp]')-options.constants(1)*eye(2*nx+npn+nz);
    prob               = iqclmi(prob,Vtot,-1,ConstMain,OuterMain);

    Vcoupling1         = blkdiag(X,-Q,Yhat,-X,oblkdiag(Be));
    OuterCoupling1     = [blkdiag([Bw;eye(nw)],[eye(nx);eye(nx)]);Bw,eye(nx);Dyw,zeros(ny,nx)];
    ConstCoupling1     = -options.constants(2)*eye(nx+nw);
    prob               = iqclmi(prob,Vcoupling1,-1,ConstCoupling1,OuterCoupling1);

    nxpsi1             = size(prob.Psi1.a,1);
    nxpsi2             = size(prob.Psi2.a,1);
    nxPl               = size(Pl.a,1);
    Yhat11             = Yhat(1:nxpsi1,1:nxpsi1);
    Yhat13             = Yhat(1:nxpsi1,nxpsi1+nxpsi2+1:nx);
    Yhat33             = Yhat(nxpsi1+nxpsi2+1:nx,nxpsi1+nxpsi2+1:nx);
    Vcoupling2         = blkdiag([Yhat11,Yhat13;Yhat13',Yhat33],-Y);
    OuterCoupling2     = [eye(nxpsi1+nxPl);eye(nxpsi1),zeros(nxpsi1,nxPl)];
    ConstCoupling2     = options.constants(3)*eye(nxpsi1+nxPl);
    prob               = iqclmi(prob,Vcoupling2,1,ConstCoupling2,OuterCoupling2);

    XQ                 = blkdiag(diag(diag(Q)),-gamma);
    E                  = ones(nw+1,1);
    prob               = iqclmi(prob,XQ,-1,0,E);

    if strcmp(options.StrProp,'no')
        Xnorm          = oblkdiag(De);
        OuterNorm      = blkdiag(eye(nz),Dyw);
        ConstNorm      = -1e-8*eye(nz+nw)-fOblkdiag(Dzw);
        prob           = iqclmi(prob,Xnorm,1,ConstNorm,OuterNorm);
    end
    

    disp('-------------------------------------------------------------------')
    disp('Solve main synthesis LMIs by minimizing gamma');
    prob               = iqcsolve(prob,gamma);
    
    if prob.gamma > 0
        ga             = sqrt(prob.gamma);
    elseif prob.gamma == -1
        ga             = -1;
    end
end

if ga == -1
    error('The LMIs are not feasible.');
end

% Compute suboptimal solution
if options.subopt > 1
    ga                 = options.subopt * ga;

    % Define LMI problem
    probOpt.FeasbRad   = options.FeasbRad;
    probOpt.Parser     = options.Parser;
    probOpt.Solver     = options.Solver;
    probOpt.Terminate  = options.Terminate;
    probOpt.RelAcc     = options.RelAcc;
    probOpt.eps        = 1e-9;
    prob               = iqcprob(probOpt);

    % Construct IQC-multipliers
    for i = 1:length(Delta)
        switch class(Delta{i})
            case 'ultis'
                PD = get(Delta{i},'PrimalDual');
                if strcmp(PD,'Primal')
                    set(Delta{i},'SampleTime',Ts);
                    prob = iqcltis_p(Delta{i},prob);
                elseif strcmp(PD,'Dual')
                    error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
                end            
            case 'ultid'
                PD = get(Delta{i},'PrimalDual');
                if strcmp(PD,'Primal')
                        set(Delta{i},'SampleTime',Ts);
                        prob = iqcltid_p(Delta{i},prob);
                elseif strcmp(PD,'Dual')
                    error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
                end
            case 'ultv'
                PD = get(Delta{i},'PrimalDual');
                if strcmp(PD,'Primal')
                    set(Delta{i},'SampleTime',Ts);
                    prob = iqcltv_p(Delta{i},prob);
                elseif strcmp(PD,'Dual')
                    error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
                end
            case 'ultv_rb'
                PD = get(Delta{i},'PrimalDual');
                if strcmp(PD,'Primal')
                    set(Delta{i},'SampleTime',Ts);
                    prob = iqcltv_rb_p(Delta{i},prob);
                elseif strcmp(PD,'Dual')
                    error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
                end
            case 'udel'
                PD = get(Delta{i},'PrimalDual');
                if strcmp(PD,'Primal')
                    Mu   = To'*M*Ti';
                    Mu   = Mu(Delta{i}.OutputChannel{1},Delta{i}.InputChannel{1}).d;
                    if norm(Mu) < 1e-14
                        set(Delta{i},'MstrictlyProp','yes');
                    else
                        set(Delta{i},'MstrictlyProp','no');
                    end
                    set(Delta{i},'SampleTime',Ts);
                    prob = iqcdel_p(Delta{i},prob);
                elseif strcmp(PD,'Dual')
                    error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
                end
            case 'usbsr'
                PD = get(Delta{i},'PrimalDual');
                if strcmp(PD,'Primal')
                    set(Delta{i},'SampleTime',Ts);
                    prob = iqcsbsr_p(Delta{i},prob);
                elseif strcmp(PD,'Dual')
                    error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
                end
            case 'usg'
                PD = get(Delta{i},'PrimalDual');
                if strcmp(PD,'Primal')
                    set(Delta{i},'SampleTime',Ts);
                    prob = iqcsg_p(Delta{i},prob);
                elseif strcmp(PD,'Dual')
                    error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
                end
            case 'ups'
                PD = get(Delta{i},'PrimalDual');
                if strcmp(PD,'Primal')
                    set(Delta{i},'SampleTime',Ts);
                    prob = iqcps_p(Delta{i},prob);
                elseif strcmp(PD,'Dual')
                    error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
                end
        end
    end

    % Construct extended plant compatible with the structure of "Pi"
    Psi                = get(prob,'Psi');
    sc                 = get(prob,'sc');
    IO                 = get(prob,'IO');
    
    Psi1               = sc*get(prob,'Psi1');
    Psi1m              = fss2m(Psi1);
    P                  = get(prob,'P');
    Y                  = iqcvar(prob,[size(Psi1.a,1),size(Psi1.a,1)],'symmetric');
    kyp                = [0,1;1,0];
    V1                 = blkdiag(kron(kyp,Y),P);
    prob               = iqclmi(prob,V1,1,options.constants(4)*eye(size(Psi1m,2)),Psi1m);
    
    % Construct extended plant compatible with the structure of "Pi"
    Plex               = fplantex(Pl,IO);
    
    Plex               = ss(Plex.a,Plex.b,[Plex.c(1:end-ny,:);zeros(nw,size(Plex.a,1));Plex.c(end-ny+1:end,:)],...
                            [Plex.d(1:end-ny,:);zeros(nw,size(Plex.d,2)-nw),eye(nw);Plex.d(end-ny+1:end,:)]);

    Psi_sc             = ss(sc*Psi);Psi_sc.Ts = Ts;
    Psiaug             = fAugss(Psi_sc,ss([],[],[],eye(nz+ny+nw),Ts),1);
    G                  = fmultss(Psiaug,Plex);
    [nGo,nGi]          = size(G);
    nx                 = size(G.a,1);

    A                  = G.a;
    B                  = G.b;
    Bp                 = G.b(:,1:end-nw);
    Bw                 = G.b(:,end-nw+1:end);
    npn                = size(Bp,2);
    
    CPsi               = G.c(1:nGo-nz-nw-ny,:);
    
    DPsi               = G.d(1:nGo-nz-nw-ny,:);
    DPsip              = G.d(1:nGo-nz-nw-ny,1:end-nw);
    DPsiw              = G.d(1:nGo-nz-nw-ny,end-nw+1:end);
    
    Cz                 = G.c(nGo-nz-nw-ny+1:1:nGo-nw-ny,:);
    Dz                 = G.d(nGo-nz-nw-ny+1:1:nGo-nw-ny,:);
    Dzp                = G.d(nGo-nz-nw-ny+1:1:nGo-nw-ny,1:end-nw);
    Dzw                = G.d(nGo-nz-nw-ny+1:1:nGo-nw-ny,end-nw+1:end);
    
    Dw                 = G.d(nGo-nw-ny+1:1:nGo-ny,:);
    Cy                 = G.c(nGo-ny+1:1:end,:);
    Dy                 = G.d(nGo-ny+1:1:end,:);
    Dyp                = G.d(nGo-ny+1:1:end,1:end-nw);
    Dyw                = G.d(nGo-ny+1:1:end,end-nw+1:end);
    
    % Define LMI variables
    Yhat               = iqcvar(prob,[nx,nx],'symmetric');
    X                  = iqcvar(prob,[nx,nx],'symmetric');
    Ae                 = iqcvar(prob,[nx,nx],'full');
    Be                 = iqcvar(prob,[nx,ny],'full');
    Ce                 = iqcvar(prob,[nz,nx],'full');
    if strcmp(options.StrProp,'no')
        De             = iqcvar(prob,[nz,ny],'full');
    elseif strcmp(options.StrProp,'yes')
        De             = zeros(nz,ny);
    end
    gamma              = iqcvar(prob,[1,1],'symmetric');

    if strcmp(options.perf,'L2')
        Xtot           = [blkdiag(Yhat,[X,Ae,Be,zeros(nx,ny)]);zeros(nGi,3*nx+2*ny);zeros(nz,2*nx),Ce,De,De];
        Vtot           = blkdiag(oblkdiag(Xtot),P);
        OuterMain      = [eye(2*nx+nGi+nz);A,A,B,zeros(nx,nz);zeros(nx),A,B,zeros(nx,nz);...
                          eye(nx),zeros(nx,nx+nGi+nz);zeros(ny,nx),Cy,Dy,zeros(ny,nz);...
                          Cy,zeros(ny,nx+nGi+nz);CPsi,CPsi,DPsi,zeros(size(CPsi,1),nz)];
        Const_intm     = [zeros(nw,2*nx),Dw,zeros(nw,nz);zeros(nz,2*nx+nGi),eye(nz)];
        ConstMain      = ga*(Const_intm'*Const_intm) - fOblkdiag([Cz,Cz,Dz]')-options.constants(1)*eye(2*nx+nGi+nz);
        prob           = iqclmi(prob,Vtot,-1,ConstMain,OuterMain);
            
        Vcoupling1     = blkdiag(X,-Yhat);
        OuterCoupling1 = [eye(nx);eye(nx)];
        ConstCoupling1 = options.constants(2)*eye(nx);
        prob           = iqclmi(prob,Vcoupling1,1,ConstCoupling1,OuterCoupling1);
        
        nxpsi1         = size(prob.Psi1.a,1);
        nxpsi2         = size(prob.Psi2.a,1);
        nxPl           = size(Pl.a,1);
        
        Yhat11         = Yhat(1:nxpsi1,1:nxpsi1);
        Yhat13         = Yhat(1:nxpsi1,nxpsi1+nxpsi2+1:nx);
        Yhat33         = Yhat(nxpsi1+nxpsi2+1:nx,nxpsi1+nxpsi2+1:nx);
        Vcoupling2     = blkdiag([Yhat11,Yhat13;Yhat13',Yhat33],-Y);
        OuterCoupling2 = [eye(nxpsi1+nxPl);eye(nxpsi1),zeros(nxpsi1,nxPl)];
        ConstCoupling2 = options.constants(3)*eye(nxpsi1+nxPl);
        prob           = iqclmi(prob,Vcoupling2,1,ConstCoupling2,OuterCoupling2);
        
        Vnorm2         = [gamma*eye(nx+nz),[Ae,Be;Ce,De];[Ae,Be;Ce,De]',gamma*eye(nx+ny)];
        prob           = iqclmi(prob,Vnorm2,1);
    
        disp('-------------------------------------------------------------------')
        disp('Minimize the matrix norm of the estimator system matrices          ');
        disp('|| [AE,BE;CE,DE] || for a suboptimal gamma. ');
    
        prob           = iqcsolve(prob,gamma);

    elseif strcmp(options.perf,'H2')
        Q              = iqcvar(prob,[nw,nw],'symmetric');
    
        Xtot           = [blkdiag(Yhat,[X,Ae,Be,zeros(nx,ny)]);zeros(npn,3*nx+2*ny);zeros(nz,2*nx),Ce,De,De];
        Vtot           = blkdiag(oblkdiag(Xtot),P);
        OuterMain      = [eye(2*nx+npn+nz);A,A,Bp,zeros(nx,nz);zeros(nx),A,Bp,zeros(nx,nz);...
                             eye(nx),zeros(nx,nx+npn+nz);zeros(ny,nx),Cy,Dyp,zeros(ny,nz);...
                             Cy,zeros(ny,nx+npn+nz);CPsi,CPsi,DPsip,zeros(size(CPsi,1),nz)];
        ConstMain      = blkdiag(zeros(2*nx+npn),eye(nz))+fOblkdiag([-Cz,-Cz,-Dzp]')-options.constants(1)*eye(2*nx+npn+nz);
        prob           = iqclmi(prob,Vtot,-1,ConstMain,OuterMain);
    
        Vcoupling1     = blkdiag(X,-Q,Yhat,-X,oblkdiag(Be));
        OuterCoupling1 = [blkdiag([Bw;eye(nw)],[eye(nx);eye(nx)]);Bw,eye(nx);Dyw,zeros(ny,nx)];
        ConstCoupling1 = -options.constants(2)*eye(nx+nw);
        prob           = iqclmi(prob,Vcoupling1,-1,ConstCoupling1,OuterCoupling1);
    
        nxpsi1         = size(prob.Psi1.a,1);
        nxpsi2         = size(prob.Psi2.a,1);
        nxPl           = size(Pl.a,1);
        Yhat11         = Yhat(1:nxpsi1,1:nxpsi1);
        Yhat13         = Yhat(1:nxpsi1,nxpsi1+nxpsi2+1:nx);
        Yhat33         = Yhat(nxpsi1+nxpsi2+1:nx,nxpsi1+nxpsi2+1:nx);
        Vcoupling2     = blkdiag([Yhat11,Yhat13;Yhat13',Yhat33],-Y);
        OuterCoupling2 = [eye(nxpsi1+nxPl);eye(nxpsi1),zeros(nxpsi1,nxPl)];
        ConstCoupling2 = options.constants(3)*eye(nxpsi1+nxPl);
        prob           = iqclmi(prob,Vcoupling2,1,ConstCoupling2,OuterCoupling2);
    
        XQ             = diag(diag(Q));
        E              = ones(nw,1);
        prob           = iqclmi(prob,XQ,-1,ga^2,E);
        
        if strcmp(options.StrProp,'no')
            Xnorm      = oblkdiag(De);
            OuterNorm  = blkdiag(eye(nz),Dyw);
            ConstNorm  = -1e-8*eye(nz+nw)-fOblkdiag(Dzw);
            prob       = iqclmi(prob,Xnorm,1,ConstNorm,OuterNorm);
        end
        
        Vnorm2         = [gamma*eye(nx+nz),[Ae,Be;Ce,De];[Ae,Be;Ce,De]',gamma*eye(nx+ny)];
        prob           = iqclmi(prob,Vnorm2,1);
    
        disp('-------------------------------------------------------------------')
        disp('Minimize the matrix norm of the estimator system matrices          ');
        disp('|| [AE,BE;CE,DE] || for a suboptimal gamma. ');
        prob           = iqcsolve(prob,gamma);
    end
end

if prob.gamma == -1
    error('The LMIs are not feasible.');
end

Yhatn                  = iqcdec2mat(prob,Yhat);
Xn                     = iqcdec2mat(prob,X);
Aen                    = iqcdec2mat(prob,Ae);
Ben                    = iqcdec2mat(prob,Be);
Cen                    = iqcdec2mat(prob,Ce);

if strcmp(options.StrProp,'no')
    Den                = iqcdec2mat(prob,De);
elseif strcmp(options.StrProp,'yes')
    Den                = zeros(nz,ny);
end

if (strcmp(options.perf,'H2')) && (norm(Den) < 1e-6)
    Den                = 0*Den;
end

AE                     = (Aen-Xn*A-Ben*Cy)/(Yhatn-Xn);
BE                     = Ben;
CE                     = Cen/(Yhatn-Xn);
DE                     = Den;
E                      = balreal(ss(AE,BE,CE,DE));
end