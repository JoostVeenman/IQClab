function prob = iqcinvariance(M,Delta,varargin)
% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        26-04-2022
% 
% -------------------------------------------------------------------------
%
% Description: This function performs IQC-based invariance analyses
%
% Syntax:      prob = iqcinvariance(M,Delta,varargin)
%
% Usage:       prob = iqcinvariance(M,Delta) performs an IQC-based
%              invariance analysis for the uncertain plant lft(Delta,M).
%              Here: 
%                1.) M is an LTI plant (continuous or discrete time)
%                2.) Delta = {Delta1,...,DeltaN,Perf1,...PerfM} is the
%                    uncertainty/performance block whose entries are
%                    objects from the class iqcdelta and which have been
%                    associated with an IQC multiplier using "iqcassign".
%                    The uncertainty/performance block Delta/P should be
%                    provided as a cell.
%
%              It is possible to specify several properties:
%
%                  prob = iqcinvariance(M,Delta,...
%                               'prop1','value1','prop2','value2',...)
%
%              The properties that can be specified are discussed in the 
%              class "iqcprob".
%
%              The algorithm currently covers 4 types of analyses:
%
%                1.) Given the uncertain plant
%
%                      q = M11 p + M12 w,   p = \Delta(q)
%
%                    with M admitting the realization
%
%                      xdot =  A x +  B1 p +  B2 w,     x(0) = 0
%                         q = C1 x + D11 p + D12 w
%
%                    and the fixed constant alp>0 (default is 1) compute
%                    (minimize) the ellipsoidal region H such that
%
%                     x(t)\in{x\in R^n: x^T(t)H x(t)\leq
%                              alp^2\|w\|_2^2} for all t\geq0           (1)
%
%                    This option can be selected by specifying the
%                    performance metric option 'e2x' (energy to state) to
%                    the performance block.
%
%                    Outputs:
%                     - H:   The ellipsoid H satisfying (1) for a given
%                            value of alp
%                     - gamma: Trace(H^-1)
%
%                2.) Given the uncertain plant
%
%                      q = M11 p + M12 w,   p = \Delta(q)
%                      z = M21 p + M22 w
%
%                    with M admitting the realization
%
%                      xdot =  A x +  B1 p +  B2 w,     x(0) = 0
%                         q = C1 x + D11 p + D12 w
%                         z = C2 x
%
%                    with M21(\infty) = M22(\infty) = 0 (i.e. M21 and M22
%                    are strictly proper) and the fixed constant alp>0
%                    (default is 1), compute (minimize) the ellipsoidal
%                    region H such that
%
%                     z(t)\in{z\in R^n_z: z^T(t)H z(t)\leq
%                              alp^2\|w\|_2^2} for all t\geq0           (2)
%
%                    This option can be selected by specifying the
%                    performance metric option 'e2z' (energy to output) to
%                    the performance block.
%
%                    Outputs:
%                     - H:   The ellipsoid H satisfying (2) for a given
%                            value of alp.
%                     - gamma: Trace(H^-1)
%
%                3.) Given the uncertain plant
%
%                      q = M11 p + M12 w,   p = \Delta(q)
%                      z = M21 p + M22 w
%
%                    with M admitting the realization
%
%                      xdot =  A x +  B1 p +  B2 w,     x(0) = 0
%                         q = C1 x + D11 p + D12 w
%                         z = C2 x
%
%                    with M21(\infty) = M22(\infty) = 0 (i.e. M21 and M22
%                    are strictly proper and the fixed constant alp>0
%                    (default is 1), compute the peak gains gam_j j \in
%                    {1,...,nz} on the performance output channels such
%                    that
%
%                    |z_j(t)| \leq sqrt(gam_j)*alp\|d\|_2 for all t\geq 0, 
%                      j\in {1,...,n}                                   (3)
%
%                    by minimizing the sum of gam_j, j \in {1,...,n}.
%
%                    This option can be selected by specifying the
%                    performance metric option e2p to the performance
%                    block.
%
%                    Output:
%                     - PeakGains = sqrt(gam_j)*alp
%
%                      (If not specified with iqcdelta it is assumed that
%                      alp = 1) 
%
%                4.) Given the uncertain plant
%
%                      q = M11 p,   p = \Delta(q)
%                      z = M21 p
%
%                    with M admitting the realization
%
%                      xdot =  A x +  B1 p,     x(0) = x0
%                         q = C1 x + D11 p
%                         z = C2 x
%
%                    with M21(\infty) = 0 (i.e. M21 is strictly proper, and
%                    the non-zero initial state x0 minimize gamma such that
%
%                      z(t)\in{z\in R^n_z: z^THz\leq alp||x0|| 
%                                                   for all t\geq0      (4)
%
%                    This option can be selected by specifying the
%                    performance metric option x2z (state to peak) to the
%                    performance block.
%
%                    Output:
%                     - H:   The ellipsoid H satisfying (4) for a given
%                            value of alp. 
%                     - gamma: Trace(H^-1)
% -------------------------------------------------------------------------

if isempty(varargin)
    prob         = iqcinvar;
elseif length(varargin) == 1
    prob         = iqcinvar;
    if isstruct(varargin{1})
        f        = fieldnames(varargin{1});
        for i = 1:length(f)
            set(prob,f{i},eval(['varargin{1}.',f{i}]))
        end
    else
        error('Error: The input arguments should come in pairs, or should be defined as a structure.');
    end
elseif length(varargin) > 1
    prob         = iqcinvar;
    if mod(length(varargin),2) ~= 0
        error('Error: The input arguments should come in pairs, or should be defined as a structure.');
    end
    for i = 1:length(varargin)/2
        set(prob,varargin{2*i-1},varargin{2*i})
    end
end

DE               = exist('Delta','var');
if DE
    if isobject(Delta)
       U{1}      = Delta;
       Delta     = U;
    end
else
    Delta        = {};
end

% Sort performance and uncertainty blocks
k                = 1;
j                = 1;
for i = 1:length(Delta)
    if isa(Delta{i},'iqcdelta')
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

DU               = exist('DeltaU','var');
DP               = exist('DeltaP','var');
if DU && DP
    Delta        = [DeltaU,DeltaP];
elseif DU && ~DP
    Delta        = DeltaU;
elseif ~DU && DP
    Delta        = DeltaP;
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
if ~isempty(To) && ~isempty(Ti)
    M             = To*M*Ti;
end

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
pm               = [0,0,0,0,0];
chi              = [];
cho              = [];

for i = 1:length(Delta)
    switch class(Delta{i})
        case 'ultis'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                set(Delta{i},'TerminalCost','on');
                prob = iqcltis_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end            
        case 'ultid'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                set(Delta{i},'TerminalCost','on');
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
                set(Delta{i},'TerminalCost','on');
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
                set(Delta{i},'TerminalCost','on');
                prob = iqcdel_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'usbsr'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                set(Delta{i},'TerminalCost','on');
                prob = iqcsbsr_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'usg'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                set(Delta{i},'TerminalCost','on');
                prob = iqcsg_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'ups'
            PD = get(Delta{i},'PrimalDual');
            if strcmp(PD,'Primal')
                set(Delta{i},'SampleTime',Ts);
                set(Delta{i},'TerminalCost','on');
                prob = iqcps_p(Delta{i},prob);
            elseif strcmp(PD,'Dual')
                error('Error: Dual IQC-multipliers are not supported for the function iqcanalysis.');
            end
        case 'iqcdelta'
            if strcmp(Delta{i}.PerfMetric,'e2x')
                pm(k,1:4) = [1,0,0,0];
            elseif strcmp(Delta{i}.PerfMetric,'e2p')
                pm(k,1:4) = [0,1,0,0];
            elseif strcmp(Delta{i}.PerfMetric,'x2z') 
                pm(k,1:4) = [0,0,1,0];
            elseif strcmp(Delta{i}.PerfMetric,'e2z') 
                pm(k,1:4) = [0,0,0,1];
            else
                error('Please consider the function iqcanalysis for this type of performance metrics');
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
if      spm(1) > 0 && spm(2) == 0 && spm(3) == 0 && spm(4) == 0
    pm  = 'e2x';
elseif spm(1) == 0 && spm(2) > 0  && spm(3) == 0 && spm(4) == 0
    pm  = 'e2p';
elseif spm(1) == 0 && spm(2) == 0 && spm(3) > 0  && spm(4) == 0
    pm  = 'x2z';
elseif spm(1) == 0 && spm(2) == 0 && spm(3) == 0 && spm(4) > 0
    pm  = 'e2z';
elseif spm(1) == 0 && spm(2) == 0 && spm(3) == 0 && spm(4) == 0
    error('Please consider the function iqcanalysis for checking robust stability.');
else
    error('Error: It is only possible to consider one single performance metric in the IQC-analyis');
end
% update performance channels
[m1,m2] = size(M);
cho     = m1 + 1 - length(cho):m1;
chi     = m2 + 1 - length(chi):m2;

% Construct extended plant compatible with the structure of "Pi"
sc      = get(prob,'sc');
Psi     = sc*get(prob,'Psi');Psi.Ts = Ts;
IO      = get(prob,'IO');
P       = get(prob,'P');

% Define KYP certificate structure (for continuous or discrete time analysis)
if Ts == 0
    kyp = [0,1;1,0];
else
    kyp = [-1,0;0,1];
end

% Define terminal cost constraint 
Z       = get(prob,'Z');
Tz      = get(prob,'Tz');

switch pm
    case 'e2x'
        % Extract alp from DeltaP 
        alp            = prob.alp;
        
        % Construct extended plant compatible with the structure of "Pi"
        Mex            = fplantex(M,IO);
        
        % Construct the main LMI matrices
        G              = fmultss(Psi,Mex);
        nx             = size(G.a,1);
        nxp            = size(Mex.a,1);
        nxpsi          = size(Psi.a,1);
        B              = fss2m(G);
        sB             = size(B,2)-length(chi);
        A              = alp^2*fJt(length(chi),sB)*fJt(length(chi),sB)';
        B1             = [eye(nx+nxp);zeros(nxpsi,nxp),Tz,zeros(nxpsi,nxp)];
        A1             = -fOblkdiag(eye(nxp),zeros(nxpsi));
        Tg             = ones(nxp+1,1);

        % Define LMI variables
        gamma          = iqcvar(prob,[1,1],'symmetric');
        X              = iqcvar(prob,[nx,nx],'symmetric');
        H              = iqcvar(prob,[nxp,nxp],'symmetric');
        dH             = diag(diag(H));
        V              = blkdiag(kron(kyp,X),P);
        V1             = blkdiag(H,X,-Z);
        Vg             = blkdiag(dH,-gamma);

        % Define LMIs
        prob           = iqclmi(prob,gamma,1);
        prob           = iqclmi(prob,V,-1,A,B);
        prob           = iqclmi(prob,V1,1,A1,B1);
        prob           = iqclmi(prob,Vg,-1,0,Tg);
        prob           = iqclmi(prob,gamma,-1,get(prob,'gmax'));
        prob           = iqcsolve(prob,gamma);
        
        if prob.gamma ~= -1
            % Check solution
            disp('The solver thinks the problem is feasible.')
            if strcmp(prob.SolChk,'on')
                eps    = 1e-14;
                
                Vn     = iqcdec2mat(prob,V);
                eP     = eig(B'*Vn*B+eps*eye(size(B,2))-A);
                
                V1n    = iqcdec2mat(prob,V1);
                eP1    = eig(B1'*V1n*B1-A1-eps*eye(size(B1,2)));
                
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0 || sum(eP1 < 0) > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                    prob.H   = -1;
                else
                    disp(['The problem is feasible: gamma = ',num2str(prob.gamma)]);
                    prob.H    = iqcdec2mat(prob,H)^-1;
                end
            else
                prob.H   = iqcdec2mat(prob,H)^-1;
            end
        else
            disp('The solver thinks the problem is infeasible.')
            prob.H   = -1;
        end
    case 'e2p'
        if norm(M.d(cho,:)) > 1e-15
            error('The uncertain plant must be strictly proper from all inputs to all performnace outputs.');
        else
            M.d(cho,:) = 0*M.d(cho,:);
        end
        % Extract alp from DeltaP 
        alp            = prob.alp;
        
        % Construct extended plant compatible with the structure of "Pi"
        Mex            = fplantex(M(1:end-length(cho),:),IO);
        
        % Construct the main LMI matrices
        G              = fmultss(Psi,Mex);
        nx             = size(G.a,1);
        nxp            = size(Mex.a,1);
        nxpsi          = size(Psi.a,1);
        B              = fss2m(G);
        sB             = size(B,2)-length(chi);
        A              = alp^2*fJt(length(chi),sB)*fJt(length(chi),sB)';
        B1             = [eye(1+nxpsi+nxp);zeros(nxpsi,1),Tz,zeros(nxpsi,nxp)];
        for i = 1:length(cho)
            A1{i}      = -fOblkdiag(M.c(cho(i),:),zeros(nxpsi));    
        end
        Tg             = ones(length(cho)+1,1);

        % Define LMI variables
        gamma          = iqcvar(prob,[1,1],'symmetric');
        X              = iqcvar(prob,[nx,nx],'symmetric');
        gn             = iqcvar(prob,[length(cho),1],'full');
        dgn            = diag(gn);
        V              = blkdiag(kron(kyp,X),P);
        for i = 1:length(cho)
            V1{i}      = blkdiag(dgn(i,i),X,-Z);
        end
        Vg             = blkdiag(dgn,-gamma);

        % Define LMIs
        prob           = iqclmi(prob,gamma,1);
        prob           = iqclmi(prob,V,-1,A,B);
        for i = 1:length(cho)
            prob       = iqclmi(prob,V1{i},1,A1{i},B1);
        end
        prob           = iqclmi(prob,Vg,-1,0,Tg);
        prob           = iqclmi(prob,gamma,-1,get(prob,'gmax'));
        prob           = iqcsolve(prob,gamma);
        
        if prob.gamma ~= -1
            % Check solution
            disp('The solver thinks the problem is feasible.')
            if strcmp(prob.SolChk,'on')
                eps    = 1e-14;
                Vn     = iqcdec2mat(prob,V);
                eP     = eig(B'*Vn*B+eps*eye(size(B,2))-A);
                eP1    = [];
                for i = 1:length(cho)
                    V1n{i} = iqcdec2mat(prob,V1{i});
                    eP1    = [eP1;eig(B1'*V1n{i}*B1-A1{i}-eps*eye(size(B1,2)))];
                end
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0 || sum(eP1 < 0) > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                    prob.PeakGain = -ones(size(gn.svar,1),1);
                else
                    disp(['The problem is feasible: gamma = ',num2str(prob.gamma)]);
                    prob.PeakGain = sqrt(iqcdec2mat(prob,gn))*alp;
                end
            else
                prob.PeakGain = sqrt(iqcdec2mat(prob,gn))*alp;
            end
        else
            disp('The solver thinks the problem is infeasible.')
            prob.PeakGain = -ones(size(gn.svar,1),1);
        end
    case 'x2z'
        if norm(M.d(cho,:)) > 1e-15
            error('The uncertain plant must be strictly proper from all inputs to all performnace outputs.');
        else
            M.d(cho,:) = 0*M.d(cho,:);
        end
        % Extract alp from DeltaP 
        alp            = prob.alp;

        % Construct extended plant compatible with the structure of "Pi"
        Mex            = fplantex(M(1:end-length(cho),:),IO);
        G              = fmultss(Psi,Mex);
        
        % Construct the main LMI matrices
        nx             = size(G.a,1);
        nxp            = size(Mex.a,1);
        nxpsi          = size(Psi.a,1);
        B              = fss2m(G);
        B1             = [eye(length(cho)+nx);zeros(nxpsi,length(cho)),Tz,zeros(nxpsi,nxp)];
        A1             = -fOblkdiag(M.c(cho,:),zeros(nxpsi));
        Tg             = ones(length(cho)+1,1);
        
        % Define LMI variables
        gamma          = iqcvar(prob,[1,1],'symmetric');
        X              = iqcvar(prob,[nx,nx],'symmetric');
        X22            = X(nxpsi+1:X.Dim(1),nxpsi+1:X.Dim(1));
        H              = iqcvar(prob,[length(cho),length(cho)],'symmetric');
        dH             = diag(diag(H));
        V              = blkdiag(kron(kyp,X),P);
        V1             = blkdiag(H,X,-Z);
        Vg             = blkdiag(dH,-gamma);

        % Define LMIs
        prob           = iqclmi(prob,gamma,1);
        prob           = iqclmi(prob,V,-1,0,B);
        prob           = iqclmi(prob,V1,1,A1,B1);
        prob           = iqclmi(prob,X22,-1,alp*eye(size(M.a,1)));
        prob           = iqclmi(prob,Vg,-1,0,Tg);
        prob           = iqclmi(prob,gamma,-1,get(prob,'gmax'));
        
        prob           = iqcsolve(prob,gamma);
        
        if prob.gamma ~= -1
            % Check solution
            disp('The solver thinks the problem is feasible.')
            if strcmp(prob.SolChk,'on')
                eps    = 1e-14;
                
                Vn     = iqcdec2mat(prob,V);
                eP     = eig(B'*Vn*B+eps*eye(size(B,2)));
                
                V1n    = iqcdec2mat(prob,V1);
                eP1    = eig(B1'*V1n*B1-A1-eps*eye(size(B1,2)));
                
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0 || sum(eP1 < 0) > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                    prob.H   = -1;
                else
                    disp(['The problem is feasible: gamma = ',num2str(prob.gamma)]);
                    prob.H   = iqcdec2mat(prob,H)^-1;
                end
            else
                prob.H   = iqcdec2mat(prob,H)^-1;
            end
        else
            disp('The solver thinks the problem is infeasible.')
            prob.H   = -1;
        end
    case 'e2z'
        if norm(M.d(cho,:)) > 1e-15
            error('The uncertain plant must be strictly proper from all inputs to all performnace outputs.');
        else
            M.d(cho,:) = 0*M.d(cho,:);
        end

        % Extract alp from DeltaP 
        alp            = prob.alp;
        
        % Construct extended plant compatible with the structure of "Pi"
        Mex            = fplantex(M(1:end-length(cho),:),IO);
        
        % Construct the main LMI matrices
        G              = fmultss(Psi,Mex);
        nx             = size(G.a,1);
        nxp            = size(Mex.a,1);
        nxpsi          = size(Psi.a,1);
        B              = fss2m(G);
        sB             = size(B,2)-length(chi);
        A              = alp^2*fJt(length(chi),sB)*fJt(length(chi),sB)';
        B1             = [eye(length(cho)+nxpsi+nxp);zeros(nxpsi,length(cho)),Tz,zeros(nxpsi,nxp)];
        A1             = -fOblkdiag(M.c(cho,:),zeros(nxpsi));
        Tg             = ones(length(cho)+1,1);

        % Define LMI variables
        gamma          = iqcvar(prob,[1,1],'symmetric');
        X              = iqcvar(prob,[nx,nx],'symmetric');
        H              = iqcvar(prob,[length(cho),length(cho)],'symmetric');
        dH             = diag(diag(H));
        V              = blkdiag(kron(kyp,X),P);
        V1             = blkdiag(H,X,-Z);
        Vg             = blkdiag(dH,-gamma);

        % Define LMIs
        prob           = iqclmi(prob,gamma,1);
        prob           = iqclmi(prob,V,-1,A,B);
        prob           = iqclmi(prob,V1,1,A1,B1);
        prob           = iqclmi(prob,Vg,-1,0,Tg);
        prob           = iqclmi(prob,gamma,-1,get(prob,'gmax'));
        prob           = iqcsolve(prob,gamma);

        if prob.gamma ~= -1
            % Check solution
            disp('The solver thinks the problem is feasible.')
            if strcmp(prob.SolChk,'on')
                eps    = 1e-14;
                
                Vn     = iqcdec2mat(prob,V);
                eP     = eig(B'*Vn*B+eps*eye(size(B,2))-A);

                V1n    = iqcdec2mat(prob,V1);
                eP1    = eig(B1'*V1n*B1-A1-eps*eye(size(B1,2)));
                
                disp('-------------------------------------------------------------------');
                disp('Perform a solution check: ...')
                if sum(eP > 0) > 0 || sum(eP1 < 0) > 0
                    disp('The problem is infeasible.');
                    set(prob,'gamma',-1);
                    prob.H   = -1;
                else
                    disp(['The problem is feasible: gamma = ',num2str(prob.gamma)]);
                    prob.H   = iqcdec2mat(prob,H)^-1;
                end
            else
                prob.H   = iqcdec2mat(prob,H)^-1;
            end
        else
            disp('The solver thinks the problem is infeasible.')
            prob.H   = -1;
        end 
end
end