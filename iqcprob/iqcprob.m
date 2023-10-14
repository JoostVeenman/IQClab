classdef iqcprob < matlab.mixin.SetGetExactNames

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
% Date:        24-11-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition to setup an IQC/LMI problem
%
% Syntax:      prob = iqcprob(varargin)
%
% Usage:       iqcprob defines an IQC-problem class, which supports the
%              parsers LMIlab as well as Yalmip.
%
%              For "prob = iqcprob(varargin)", the varargin inputs come in
%              pairs and can be defined as: 
%
%              prob = iqcprob('name','prop1','value1','prop2','value2',...)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              iqcprobOpt.prop1 = value1
%              iqcprobOpt.prop2 = value2
%                       ...              
%              iqcprobOpt.propN = valueN
% 
%              prob = iqcprob(iqcprobOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%                   set(prob,'propX','valueX') and get(prob,'propX')
%
%              The properties that can be specified are:
%            
%                1.) 'Parser' specifies either one of the following
%                    parsers: 'LMIlab' (default) and 'Yalmip'.
%
%                2.) 'Solver' specifies the Yalmip solver, such as e.g.
%                    'sdpt3' (default) or 'sedumi'. The reader is referred
%                    to https://yalmip.github.io/allsolvers/ for an
%                    overview of solvers that are supported.
%
%                3.) 'gmax' specifies the maximum worst-case performance
%                    value (default: 1e4). 
%
%                4.) 'RelAcc' is an LMIlab option that specifies the
%                    relative accuracy on the optimal value (default: 1e-6)
%
%                5.) 'MaxNumIter' is an LMIlab option that specifies the
%                    maximum number of iterations that can be performed by
%                    the optimization procedure (default: 800).
%
%                6.) 'FeasbRad' specifies the feasibility radius of the LMI
%                    problem (default 1e9).
%
%                7.) 'Terminate' is an LMIlab option that helps speed up
%                    the termination. If set to an integer value J > 0, the
%                    code terminates when the objective c^T x has not
%                    decreased by more than the desired relative accuracy
%                    during the last J iterations. 
%
%                8.) 'Display' allows to turn on/off the trace of execution
%                    of the optimization procedure (default: 'off').
%
%                9.) 'Init' can be used to specify an initial condition for
%                    the LMI optimization problem (default: []). 
%
%               10.) 'eps' allows to enforce strict LMI by a small
%                    nonnegative constant (default: 1e-9).
%
%               11.) 'SolChk' allows to check whether an obtained LMI
%                    solution is indeed feasible.
%
%               12.) 'Pi11pos' allows to enforce strictness of the
%                    constraint Pi11(iw) > 0 for all w\in R\cup{\infty}. If
%                    set to a value larger than 0, the IQC-analysis will
%                    include the constraint:
%
%                    Pi11(iw) > Pi11pos*I for all w\in R\cup{\infty}.
%
%              The properties that are obtained by an IQC-analysis are
%
%                1.) 'gamma' is the optimal value obtained by the LMI
%                    optimization problem. In case the specified
%                    performance metric is 'P' (passive) or when only
%                    performing a robust stability analysis, then
%                    feasibility is denoted by a 1. Infeasibility is always
%                    denoted by -1.
%
%                2.) 'xsol' is the LMIlab solution returned by mincx.
%
%                3.) 'lmi' is the set of LMI constraints.
%
%                4.) 'Psi' is basis function of the overall IQC-multiplier
%                    [sc*Psi(iw)]^* P[sc*Psi(iw)]. See method "fPsi" below
%                    for further info.
%
%                5.) 'P' is the overall matrix variable of the overall
%                    IQC-multiplier [sc*Psi(iw)]^* P[sc*Psi(iw)]. See
%                    method "fP" below for further info.
%
%                6.) 'sc' is a scalings matrix of the overall
%                    IQC-multiplier [sc*Psi(iw)]^* P[sc*Psi(iw)]. See
%                    method "fsc" below for further info.
%
%                7.) 'IO' contains the in-/output data of the
%                    IQC-multipliers (See method IO below for further info
%
% -------------------------------------------------------------------------

properties
    Parser     string {mustBeMember(Parser,{'LMIlab','Yalmip'})} = 'LMIlab'; % Optional Parsers: LMIlab, Yalmip
    Solver                                                       = 'mincx';  % Optional Yalmip solvers: sdpt3, sedumi (see "https://yalmip.github.io/solver" for more options)
    gmax       double {mustBeReal,mustBeFinite}                  = 1e4;      % specify a maximum value for gamma
    RelAcc     double {mustBeReal,mustBeFinite}                  = 1e-6;     % relative accuracy
    MaxNumIter double {mustBeReal,mustBeFinite}                  = 800;      % maximum number of iterations
    FeasbRad   double {mustBeReal,mustBeFinite}                  = 1e9;      % Maximum bound on variables
    Terminate  double {mustBeReal,mustBeFinite}                  = 0;        % integer > 0 speeds up termination (see options mincx for further details)
    Display    string {mustBeMember(Display,{'on','off'})}       = 'off';    % Display0 = 'on' or 'off'
    Init       double {mustBeReal,mustBeFinite}                  = [];       % Initial condition LMIs
    eps        double {mustBeReal,mustBeFinite}                  = 1e-9;     % epsilon constant to enforce strict LMIs in Yalmip
    SolChk     string {mustBeMember(SolChk,{'on','off'})}        = 'off';    % Check if LMIs are feasible (valid for function: iqcanalysis)
    BoundVars  double {mustBeReal,mustBeFinite}                  = 0;        % 0 = 'off', 1 = enforce bound on vars, 2 = min bounds on vars
    Pi11pos    double {mustBeReal,mustBeFinite}                  = 0;        % include the constraint Pi11(iw) >  Pi11pos*I for all w (default = 0)
end

% Internal properties
properties
    gamma                                                        = [];       % Optimal value of minimized variable gamma
    xsol                                                         = [];       % LMI solution
    lmi                                                          = [];       % Empty set of LMIs
    Psi                                                          = ss([]);   % Basis function IQC multiplier
    Psi1                                                         = ss([]);   % First column of basis function IQC multiplier
    Psi2                                                         = ss([]);   % Second column of basis function IQC multiplier
    IOpsi1                                                       = [];       % IO data of first column of basis function IQC multiplier
    IOpsi2                                                       = [];       % IO data of second column of basis function IQC multiplier
    P                                                            = [];       % Overall muliplier variable P
    sc                                                           = [];       % scale factor Delta
    IO                                                           = [];       % IOdata IQC multipliers
end
methods
    function obj  = iqcprob(varargin)
        if nargin == 1
            if ~isstruct(varargin{1})
                error('Error: The second input argument should be defined by a structure (see "help iqcprob" for further details)');
            end
            if isfield(varargin{1},'Parser')
                obj.Parser = varargin{1}.Parser;
            end
            if isfield(varargin{1},'Solver')
                if ischar(varargin{1}.Solver)
                    obj.Solver = varargin{1}.Solver;
                else
                    error('Error: The solver (Yalmip option) should be specified as a string (Default = sdpt3).');
                end
            end
            if isfield(varargin{1},'RelAcc')
                if isreal(varargin{1}.RelAcc) && isscalar(varargin{1}.RelAcc) && varargin{1}.RelAcc > 0
                    obj.RelAcc = varargin{1}.RelAcc;
                else
                    error('Error: The option "RelAcc" (LMIlab option) should be defined as a real positive (small: << 1) scalar (Default = 1e-6).');
                end
            end
            if isfield(varargin{1},'MaxNumIter')
                if isreal(varargin{1}.MaxNumIter) && isscalar(varargin{1}.MaxNumIter) && varargin{1}.MaxNumIter > 0
                    obj.MaxNumIter = varargin{1}.MaxNumIter;
                else
                    error('Error: The option "MaxNumIter" (LMIlab option) should be defined as a real positive (small: << 1) scalar.');
                end
            end
            if isfield(varargin{1},'FeasbRad')
                if isreal(varargin{1}.FeasbRad) && isscalar(varargin{1}.FeasbRad) && varargin{1}.FeasbRad > 0
                    obj.FeasbRad = varargin{1}.FeasbRad;
                else
                    error('Error: The option "FeasbRad" should be defined as a real positive scalar (Default = 1e9).');
                end
            end
            if isfield(varargin{1},'Terminate')
                if isreal(varargin{1}.Terminate) && isscalar(varargin{1}.Terminate) && varargin{1}.Terminate >= 0
                    obj.Terminate = varargin{1}.Terminate;
                else
                    error('Error: The option "Terminate" (LMIlab option) should be defined as a real positive scalar (Default = 1e9).');
                end
            end
            if isfield(varargin{1},'Display')
                obj.Display = varargin{1}.Display;
            end
            if isfield(varargin{1},'Init')
                if isreal(varargin{1}.Init) && ismatrix(varargin{1}.Init)
                    obj.Init = varargin{1}.Init;
                else
                    error('Error: The option "Init" (LMIlab option) should be defined as a real vector (Default = []).');
                end
            end
            if isfield(varargin{1},'eps')
                if isreal(varargin{1}.eps) && isscalar(varargin{1}.eps) && varargin{1}.eps >= 0
                    obj.eps = varargin{1}.eps;
                else
                    error('Error: The option "eps" should be defined as a real non-negative(small: << 1) scalar (Default = 1e-9).');
                end
            end
            if isfield(varargin{1},'gmax')
                if isreal(varargin{1}.gmax) && isscalar(varargin{1}.gmax) && varargin{1}.gmax > 0
                    obj.gmax = varargin{1}.gmax;
                else
                    error('Error: The option "gmax" should be defined as a real positive scalar (Default = 1e4).');
                end
            end
            if isfield(varargin{1},'SolChk')
                obj.SolChk = varargin{1}.SolChk;
            end
            if isfield(varargin{1},'Pi11pos')
                if isreal(varargin{1}.Pi11pos) && isscalar(varargin{1}.Pi11pos) && varargin{1}.Pi11pos > 0
                    obj.Pi11pos = varargin{1}.Pi11pos;
                else
                    error('Error: The option "Pi11pos" should be defined as a real non-negative (small: << 1) scalar (Default = 0).');
                end
            end
        elseif nargin > 1        
            n = length(varargin);
            if mod(n,2) ~= 0
                error('Error: The input arguments should come in pairs')
            else
                j = linspace(2,n,n/2);
            end
            for i = 1:length(j)
                prop = varargin{j(i)-1};
                switch prop
                    case 'Parser'
                        obj.Parser = varargin{j(i)};
                    case 'Solver'
                        obj.Solver = varargin{j(i)};
                        if ischar(varargin{j(i)})
                            obj.Solver = varargin{j(i)};
                        else
                            error('Error: The solver (Yalmip option) should be specified as a string (Default = sdpt3).');
                        end
                    case 'RelAcc'
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                            obj.RelAcc = varargin{j(i)};
                        else
                            error('Error: The option "RelAcc" (LMIlab option) should be defined as a real positive (small: << 1) scalar (Default = 1e-6).');
                        end
                    case 'MaxNumIter'
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                            obj.MaxNumIter = varargin{j(i)};
                        else
                            error('Error: The option "MaxNumIter" (LMIlab option) should be defined as a real positive (small: << 1) scalar.');
                        end
                    case 'FeasbRad'
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                            obj.FeasbRad = varargin{j(i)};
                        else
                            error('Error: The option "FeasbRad" should be defined as a real positive scalar (Default = 1e9).');
                        end
                    case 'Terminate'
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} >= 0
                            obj.Terminate = varargin{j(i)};
                        else
                            error('Error: The option "Terminate" (LMIlab option) should be defined as a real positive scalar (Default = 1e9).');
                        end
                    case 'Display'
                        obj.Display = varargin{j(i)};
                    case 'Init'
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                            obj.Init = varargin{j(i)};
                        else
                            error('Error: The option "Init" (LMIlab option) should be defined as a real vector (Default = []).');
                        end
                    case 'eps'
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} >= 0
                            obj.eps = varargin{j(i)};
                        else
                            error('Error: The option "eps" should be defined as a real non-negative (small: << 1) scalar (Default = 1e-9).');
                        end
                    case 'gmax'
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                            obj.gmax = varargin{j(i)};
                        else
                            error('Error: The option "gmax" should be defined as a real positive scalar (Default = 1e4).');
                        end
                    case 'SolChk'
                        obj.SolChk = varargin{j(i)};
                    case 'Pi11pos'
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} >= 0
                            obj.Pi11pos = varargin{j(i)};
                        else
                            error('Error: The option "Pi11pos" should be defined as a real non-negative (small: << 1) scalar (Default = 0).');
                        end
                end
            end
        end
        % Initialize the LMI problem
        clear yalmip
        obj.lmi = [];
        setlmis([]);
    end
    function iqcprob = fP(iqcprob,P)
    % ---------------------------------------------------------------------
    % Description: Augment the multiplier matrix variable P
    %
    % Syntax:      iqcprob = fP(iqcprob,P)
    %
    % Usage:       Specify the IQC problem "iqcprob" as input as well as
    %              the to-be-augmented matrix variable P.
    % ---------------------------------------------------------------------
        Pprev         = get(iqcprob,'P');
        if ~isobject(Pprev)
            if isempty(Pprev)
                Pprev = iqcvar(iqcprob,[0,0],'symmetric');
            end
        end
        set(iqcprob,'P',blkdiag(Pprev,P));
    end
    function iqcprob = fsc(iqcprob,sc)
    % ---------------------------------------------------------------------
    % Description: Augment the multiplier scaling matrix sc
    %
    % Syntax:      iqcprob = fsc(iqcprob,sc)
    %
    % Usage:       Specify the IQC problem "iqcprob" as input as well as
    %              the to-be-augmented multiplier scaling matrix sc.
    % ---------------------------------------------------------------------
        scprev       = get(iqcprob,'sc');
        set(iqcprob,'sc',blkdiag(scprev,sc));
    end
    function iqcprob = fIO(iqcprob,IO)
    % ---------------------------------------------------------------------
    % Description: Augment the IO multiplier data
    %
    % Syntax:      iqcprob = fIO(iqcprob,IO)
    %
    % Usage:       Specify the IQC problem "iqcprob" as input as well as
    %              the to-be-augmented IO multiplier data IO.
    %
    %              Here we recall that the multiplier Pi is factorized as
    %
    %                               nP11   nP22  nsc1  nsc2  nuo   nui
    %                               <-->   <-->  <-->  <-->  <-->  <--> 
    %              Pi = [x]^*[x]^T[ P11  , P12 ][ T11, T12 ][Psi1, 0   ]
    %                             [ P12^T, P22 ][ T21, T22 ][0   , Psi2]
    %
    %              where the dimensional info should be specified as a row-
    %              vector:
    %
    %              IO = [nuo,nui,nuo_ltv_rb,nui_ltv_rb,... 
    %                                   ...,nsc1,nsc2,np11,np22,ns1,ns2]
    %
    %              Here:
    %               - nuo        = nr of unc. outputs
    %               - nui        = nr of unc. inputs
    %               - nuo_ltv_rb = nr of unc. outputs for "ultv_rb" class
    %               - nui_ltv_rb = nr of unc. inputs for "ultv_rb" class
    %               - nP11       = nr of row/colums of P11
    %               - nP22       = nr of row/colums of P22
    %               - nsc1       = nr of columns of T11
    %               - nsc2       = nr of columns of T22
    %               - ns1        = nr of states Psi1
    %               - ns2        = nr of states Psi2
    % ---------------------------------------------------------------------
        IO_prev      = get(iqcprob,'IO');
        set(iqcprob,'IO',[IO_prev;IO]);
    end
    function iqcprob = fPsi(iqcprob,Phi1,Phi2)
    % ---------------------------------------------------------------------
    % Description: Augment the IQC multiplier outer factor
    %
    % Syntax:      iqcprob = fPsi(iqcprob,Phi1,Phi2)
    %
    % Usage:       Specify the IQC problem "iqcprob" as input as well as
    %              the to-be-augmented basis functions Phi1 and Phi2.
    % ---------------------------------------------------------------------
    
        % Retrieve current Psi, Psi1 and Psi2
        Psi_prev     = get(iqcprob,'Psi');
        Psi1_prev    = get(iqcprob,'Psi1');
        Psi2_prev    = get(iqcprob,'Psi2');

        % Update Psi, Psi1 and Psi2
        Psi_new      = fAugss(Psi_prev,Phi1,1);
        Psi_new      = fAugss(Psi_new,Phi2,1);
        Psi1_new     = fAugss(Psi1_prev,[Phi1;zeros(size(Phi2.d,1),size(Phi1.d,2))],1);
        Psi2_new     = fAugss(Psi2_prev,[zeros(size(Phi1.d,1),size(Phi2.d,2));Phi2],1);

        % Set new Psi, Psi1 and Psi2
        set(iqcprob,'Psi',Psi_new);
        set(iqcprob,'Psi1',Psi1_new);
        set(iqcprob,'Psi2',Psi2_new);

        % Retrieve current IOpsi1 and IOpsi2
        IOpsi1_prev  = get(iqcprob,'IOpsi1');
        IOpsi2_prev  = get(iqcprob,'IOpsi2');

        % Update IOpsi1 and IOpsi2
        IOpsi1_new   = [IOpsi1_prev;size(Phi1.a,1),size(Phi1.d,1),size(Phi1.d,2)];
        IOpsi2_new   = [IOpsi2_prev;size(Phi2.a,1),size(Phi2.d,1),size(Phi2.d,2)];

        % Set new Psi, IOpsi1 and IOpsi2
        set(iqcprob,'IOpsi1',IOpsi1_new);
        set(iqcprob,'IOpsi2',IOpsi2_new);
    end
end
end