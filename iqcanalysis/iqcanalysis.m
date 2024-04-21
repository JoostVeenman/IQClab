function prob = iqcanalysis(M,Delta,varargin)
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
% Date:        23-12-2019
% 
% -------------------------------------------------------------------------
%
% Description: This function performs an IQC-analysis
%
% Syntax:      prob = iqcanalysis(M,Delta,varargin)
%
% Usage:       prob = iqcanalysis(M,Delta) performs an IQC-analysis for the
%              uncertain plant lft(Delta,M). Here:
%                1.) M is a stable LTI plant (continuous or discrete time)
%                2.) Delta = {Delta1,...,DeltaN,Perf1,...PerfM} is the
%                    uncertainty/performance block whose entries are
%                    objects from the class iqcdelta and which have been
%                    associated with an IQC multiplier using "iqcassign".
%                    The uncertainty/performance block Delta/P should be
%                    provided as a cell.
%
%              It is possible to specify several properties:
%
%                  prob = iqcanalysis(M,Delta,...
%                               'prop1','value1','prop2','value2',...)
%
%              The properties that can be specified are discussed in
%              the function "iqcprob".
%
%              The algorithms covers 5 types of analyses:
%
%                1.) Robust performance analysis with induced L2-gain as
%                    perfomrance objective.
%                2.) Robust performance analysis with H2 as performance
%                    objective. 
%                3.) Robust performance analysis with generalized H2 as
%                    performance objective. 
%                4.) Robust performance analysis with passivity as
%                    performance objective.
%                5.) Robust stability analysis.
%
% -------------------------------------------------------------------------

if isempty(varargin)
    prob         = iqcprob;
elseif length(varargin) == 1
    prob         = iqcprob(varargin{1});
elseif length(varargin) > 1
    te           = ['prob = iqcprob(varargin{',num2str(1),'}'];
    for i = 2:length(varargin)
        te       = [te,',varargin{',num2str(i),'}'];
        if i == length(varargin)
            te   = [te,');'];
        end
    end
    eval(te)
end

if isobject(Delta)
   U{1}          = Delta;
   Delta         = U;
end

% Sort performance and uncertainty blocks
k                = 1;
j                = 1;
for i = 1:length(Delta)
    if strcmp(class(Delta{i}),'iqcdelta')
        if strcmp(Delta{i}.ChannelClass{1},'P')
            DeltaP{k} = Delta{i};
            k         = k + 1;
        elseif strcmp(Delta{i}.ChannelClass{1},'C')
            error('Error: It is not permitted to have open-loop control channels in the IQC analysis problem.');
        end
    else
       DeltaU{j} = Delta{i};
       j         = j + 1;
    end
end
if k == 1
    Delta        = DeltaU;
else
    Delta        = [DeltaU,DeltaP];
end

% Permute plant in- and output channels according to the delta structure
[mo,mi]          = size(M);
To               = [];
Ti               = [];
for i = 1:length(Delta)
    [Toi,Tii]    = fT(mo,mi,Delta{i}.OutputChannel,Delta{i}.InputChannel);
    To           = [To;Toi];
    Ti           = [Ti,Tii];
end
if sum(sum(To) > 1) || sum(sum(To.') > 1) || sum(sum(Ti) > 1) || sum(sum(Ti.') > 1)
   error('Error: At least 1 plant in- or output channel was assigned twice or more.');
end
M                = To*M*Ti;

% Get sample time plant M
Ts               = M.Ts;

% Test if M is nominally stable
eM               = eig(M);
if (Ts == 0) && (sum(real(eM) >= 0)) 
    error('Error: The nominal plant is not stable.');
elseif (Ts > 0) && (sum(abs(eM) >= 1))
    error('Error: The nominal plant is not stable.');
elseif (Ts == -1) && (sum(abs(eM) >= 1))
    error('Error: The nominal plant is not stable.');
end

disp('-------------------------------------------------------------------');
disp(['Uncertain plant with ' num2str(size(M,2)) ' inputs, ' num2str(size(M,1)) ' outputs, and ', num2str(size(M.a,1)),' states.']);

% Construct IQC-multipliers
k                = 1;
pm               = [0,0,0,0];
chi              = [];
cho              = [];

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
        case 'iqcdelta'
            if strcmp(Delta{i}.PerfMetric,'L2')
                pm(k,1:4) = [1,0,0,0];
            elseif strcmp(Delta{i}.PerfMetric,'H2')
                pm(k,1:4) = [0,1,0,0];
            elseif strcmp(Delta{i}.PerfMetric,'genH2')
                pm(k,1:4) = [0,0,1,0];
            elseif strcmp(Delta{i}.PerfMetric,'P')
                pm(k,1:4) = [0,0,0,1];
            elseif strcmp(Delta{i}.PerfMetric,'e2x') || strcmp(Delta{i}.PerfMetric,'e2p') || strcmp(Delta{i}.PerfMetric,'x2p') || strcmp(Delta{i}.PerfMetric,'e2z') 
                error('Please consider the function "iqcinvariance" for this type of performance metrics');
            end
            for ijk = 1:length(Delta{i}.InputChannel)
                chi = [chi,Delta{i}.InputChannel{ijk}];
            end
            for ijk = 1:length(Delta{i}.OutputChannel)
                cho = [cho,Delta{i}.OutputChannel{ijk}];
            end
            k = k + 1;
    end
end
if size(pm,1) > 1
    spm = sum(pm);
else
    spm = pm;
end
if spm(1) > 0      && spm(2) == 0 && spm(3) == 0 && spm(4) == 0
    pm  = 'L2';
elseif spm(1) == 0 && spm(2) > 0  && spm(3) == 0 && spm(4) == 0
    pm  = 'H2';
elseif spm(1) == 0 && spm(2) == 0 && spm(3) > 0  && spm(4) == 0
    pm  = 'genH2';    
elseif spm(1) == 0 && spm(2) == 0 && spm(3) == 0 && spm(4) > 0
    pm  = 'P';
elseif spm(1) == 0 && spm(2) == 0 && spm(3) == 0 && spm(4) == 0
    pm  = 'RS';
else
    error('Error: It is only possible to consider one single performance metric in the IQC-analyis');
end
if strcmp(pm,'H2') || strcmp(pm,'genH2')
    if norm(M.d(:,end-length(chi)+1:end)) > 1e-15
        error('When considering the (generalized) H2-performance metric, the uncertain plant must be strictly proper from the generalized distrubance input to all outputs.');
    end
end

% Construct extended plant compatible with the structure of "Pi"
Psi             = get(prob,'Psi');
sc              = get(prob,'sc');
IO              = get(prob,'IO');

% Define KYP certificate structure (for continuous or discrete time analysis)
if Ts == 0
    kyp         = [0,1;1,0];
else
    kyp         = [-1,0;0,1];
end

if prob.Pi11pos > 0
    Psi1        = sc*get(prob,'Psi1');
    Psi1m       = fss2m(Psi1);
    P11         = get(prob,'P');
    X11         = iqcvar(prob,[size(Psi1.a,1),size(Psi1.a,1)],'symmetric');
    V11         = blkdiag(kron(kyp,X11),P11);
    prob        = iqclmi(prob,V11,1,prob.Pi11pos*eye(size(Psi1m,2)),Psi1m);
end

% Construct extended plant compatible with the structure of "Pi"
Mex             = fplantex(M,IO);

switch pm
    case 'L2'
        % Construct the main LMI matrices
        Psi_sc  = ss(sc*Psi);Psi_sc.Ts = Ts;
        Psiaug  = fAugss(Psi_sc,ss([],[],[],eye(length(cho)),Ts),1);
        G       = fmultss(Psiaug,Mex);
        Gm      = fss2m(G);
        sGm1    = size(Gm,1)-length(cho);
        sGm2    = size(Gm,2)-length(chi);
        nx      = size(G.a,1);
        A       = blkdiag([Gm(1:sGm1,:);fJt(length(chi),sGm2)'],eye(length(cho)));
        B       = fOblkdiag(-Gm(sGm1+1:end,:)');

        % Define LMI variables
        gamma   = iqcvar(prob,[1,1],'symmetric');
        X       = iqcvar(prob,[nx,nx],'symmetric');
        P       = get(prob,'P');
        V       = blkdiag(kron(kyp,X),P,-gamma*eye(length(chi)+length(cho)));

        % Define LMIs
        prob    = iqclmi(prob,V,-1,B,A);
        prob    = iqclmi(prob,gamma,-1,get(prob,'gmax'));
        prob    = iqcsolve(prob,gamma);

        if prob.gamma ~= -1
            % Check solution
            disp('The solver thinks the problem is feasible.')
            if strcmp(prob.SolChk,'on')
                Vn  = iqcdec2mat(prob,V);
                sA2 = size(A,2);
                eps = 1e-13;
                eP  = eig(A'*Vn*A+eps*eye(sA2)-B);
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                else
                    disp(['The problem is feasible: gamma = ',num2str(prob.gamma)]);
                end
            end
        else
            disp('The solver thinks the problem is infeasible.')
        end
    case 'H2'
        % Construct the main LMI matrices
        Psi_sc  = ss(sc*Psi);Psi_sc.Ts = Ts;
        Psiaug  = fAugss(Psi_sc,ss([],[],[],eye(length(cho)),Ts),1);
        G       = fmultss(Psiaug,Mex);
        Gm      = fss2m(G);
        sGm1    = size(Gm,1)-length(cho);
        sGm2    = size(Gm,2)-length(chi);
        nx      = size(G.a,1);
        A       = Gm(1:sGm1,1:sGm2);
        B       = -Gm(sGm1+1:end,1:sGm2).'*Gm(sGm1+1:end,1:sGm2);
        Ap      = [Gm(nx+1:2*nx,sGm2+1:end);eye(length(chi))];
        Bp      = Gm(nx+1:2*nx,sGm2+1:end);
        Aq      = ones(length(chi)+1,1);

        % Define LMI variables
        gamma   = iqcvar(prob,[1,1],'symmetric');
        X       = iqcvar(prob,[nx,nx],'symmetric');
        Y       = iqcvar(prob,[length(chi),length(chi)],'symmetric');
        Z       = diag(diag(Y)); % check!
        P       = get(prob,'P');
        V       = blkdiag(kron(kyp,X),P);
        W       = blkdiag(X,-Y);
        Q       = blkdiag(Z,-gamma);
        
        % Define LMI constraints
        prob    = iqclmi(prob,V,-1,B,A);
        prob    = iqclmi(prob,W,-1,0,Ap);
        prob    = iqclmi(prob,Q,-1,0,Aq);
        prob    = iqclmi(prob,gamma,1,0);
        prob    = iqclmi(prob,gamma,-1,get(prob,'gmax'));
        prob    = iqcsolve(prob,gamma);
        
        if prob.gamma ~= -1
            disp('The solver thinks the problem is feasible.')
            % Check solution
            if strcmp(prob.SolChk,'on')
                Vn  = iqcdec2mat(prob,V);
                Xn  = iqcdec2mat(prob,X);
                gn  = iqcdec2mat(prob,gamma);
                sA2 = size(A,2);
                eps = 1e-13;
                eP  = eig(A'*Vn*A+eps*eye(sA2)-B);
                eQ  = trace(Bp'*Xn*Bp)+eps-gn;
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0 || eQ > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                else
                    disp(['The problem is feasible: gamma = ',num2str(sqrt(prob.gamma))]);
                    set(prob,'gamma',sqrt(gn));
                end
            else
                set(prob,'gamma',sqrt(prob.gamma));
            end
        else
            disp('The solver thinks the problem is infeasible.')
        end
    case 'genH2'
        % Construct the main LMI matrices
        Psi_sc  = ss(sc*Psi);Psi_sc.Ts = Ts;
        Psiaug  = fAugss(Psi_sc,ss([],[],[],eye(length(cho)),Ts),1);
        G       = fmultss(Psiaug,Mex);
        Gm      = fss2m(G);
        sGm1    = size(Gm,1)-length(cho);
        sGm2    = size(Gm,2)-length(chi);
        nx      = size(G.a,1);
        A       = Gm(1:sGm1,1:sGm2);
        B       = -Gm(sGm1+1:end,1:sGm2).'*Gm(sGm1+1:end,1:sGm2);
        Ap      = [Gm(nx+1:2*nx,sGm2+1:end);eye(length(chi))];
        Bp      = Gm(nx+1:2*nx,sGm2+1:end);

        % Define LMI variables
        gamma   = iqcvar(prob,[1,1],'symmetric');
        X       = iqcvar(prob,[nx,nx],'symmetric');
        P       = get(prob,'P');
        V       = blkdiag(kron(kyp,X),P);
        W       = blkdiag(X,-gamma*eye(length(chi)));
        
        % Define LMI constraints
        prob    = iqclmi(prob,V,-1,B,A);
        prob    = iqclmi(prob,W,-1,0,Ap);
        prob    = iqclmi(prob,gamma,1,0);
        prob    = iqclmi(prob,gamma,-1,get(prob,'gmax'));
        prob    = iqcsolve(prob,gamma);
        
        if prob.gamma ~= -1
            disp('The solver thinks the problem is feasible.')
            % Check solution
            if strcmp(prob.SolChk,'on')
                Vn  = iqcdec2mat(prob,V);
                Xn  = iqcdec2mat(prob,X);
                gn  = iqcdec2mat(prob,gamma);
                sA2 = size(A,2);
                eps = 1e-13;
                eP  = eig(A'*Vn*A+eps*eye(sA2)-B);
                eQ  = eig(Bp'*Xn*Bp+eps-gn*eye(length(chi)));
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0 || sum(eQ > 0) > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                else
                    disp(['The problem is feasible: gamma = ',num2str(sqrt(prob.gamma))]);
                    set(prob,'gamma',sqrt(gn));
                end
            else
                set(prob,'gamma',sqrt(prob.gamma));
            end
        else
            disp('The solver thinks the problem is infeasible.')
        end
    case 'P'
        % Construct the main LMI matrices
        if length(cho) ~= length(chi)
            error('Error: The number of performance in- and outputs should be equal.');
        end
        Psi_sc  = ss(sc*Psi);Psi_sc.Ts = Ts;
        Psiaug  = fAugss(Psi_sc,ss([],[],[],eye(length(cho)),Ts),1);
        G       = fmultss(Psiaug,Mex);
        Gm      = fss2m(G);
        sGm1    = size(Gm,1)-length(cho);
        sGm2    = size(Gm,2)-length(chi);
        nx      = size(G.a,1);
        Gm1     = Gm(1:sGm1,:);
        Gm2     = [Gm(sGm1+1:end,:);fJt(length(chi),sGm2)'];
        A       = Gm1;
        B       = -Gm2.'*fOblkdiag(-eye(length(cho)))*Gm2;
        
        % Define LMI variables
        X       = iqcvar(prob,[nx,nx],'symmetric');
        P       = get(prob,'P');
        V       = blkdiag(kron(kyp,X),P);

        % Define LMI constraints
        prob    = iqclmi(prob,V,-1,B,A);
        prob    = iqcsolve(prob);
        
        if strcmp(prob.gamma,'feasible')
            disp('The solver thinks the problem is feasible.')
            % Check solution
            if strcmp(prob.SolChk,'on')
                Vn  = iqcdec2mat(prob,V);
                sA2 = size(A,2);
                eps = 1e-13;
                eP  = eig(A'*Vn*A+eps*eye(sA2)-B);
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                else
                    disp('The problem is feasible: The system is strictly (robustly) passive');
                    set(prob,'gamma',1);
                end
            end
        else
            disp('The problem is infeasible.')
        end
    case 'RS'
        % Construct the main LMI matrices
        Psi_sc  = ss(sc*Psi);Psi_sc.Ts = Ts;
        G       = fmultss(Psi_sc,Mex);
        A       = fss2m(G);
        
        % Define LMI variables
        X       = iqcvar(prob,[size(G.a,1),size(G.a,1)],'symmetric');
        P       = get(prob,'P');
        V       = blkdiag(kron(kyp,X),P);
        
        % Define LMIs
        prob    = iqclmi(prob,V,-1,0,A);
        prob    = iqcsolve(prob);
        
        if strcmp(prob.gamma,'feasible')
            disp('The solver thinks the problem is feasible.')
            % Check solution
            if strcmp(prob.SolChk,'on')
                Vn  = iqcdec2mat(prob,V);
                sA2 = size(A,2);
                eps = 1e-13;
                eP  = eig(A'*Vn*A+eps*eye(sA2));
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                else
                    disp('The problem is feasible.');
                    set(prob,'gamma',1);
                end
            end
        else
            disp('The solver thinks the problem is infeasible.')
        end
end
end