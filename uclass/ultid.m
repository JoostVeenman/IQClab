classdef ultid < handle & matlab.mixin.SetGetExactNames

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
% Date:        05-11-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition for LTI diagonally repeated dynamic
%              uncertainties of the form:
%
%                                         [ delta ,     ,  0   ]
%              Delta = Kron(I_nr,delta) = [       , ... ,      ]
%                                         [   0   ,     , delta]
%
%              where:
%              - delta \in RH_\infty
%              - delta has dimension n_out x n_in
%              - I_nr is the number of repetitions of the delta block
%              - ||delta||_\infty < alpha
%
% Syntax:      delta = ultid('name')
%              delta = ultid('name',varargin)
%
% Usage:       "delta = ultid('name')" defines a SISO LTI dynamic
%              uncertainty, which is repeated once and satisfies
%              ||delta||_\infty < 1. 
%
%              For "delta = ultid('name',varargin)", the varargin inputs
%              come in pairs and can be defined as: 
%
%              delta = ultid('name','prop1','value1','prop2','value2',...)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              ultidOpt.prop1 = value1
%              ultidOpt.prop2 = value2
%                       ...              
%              ultidOpt.propN = valueN
% 
%              delta = ultid('name',ultidOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%              set(delta,'propX','valueX') and get(delta,'propX').
%
%              The properties that can be specified are:
%            
%               1.) 'Dimensions' specifies the size of the uncertainty
%                   delta as [n_out, n_in] (default = [1,1]).
%
%               2.) 'NormBounds' specifies the H_\infty norm of the
%                   uncertainty by means of the value alpha (default = 1).
%
%               3.) 'NumberOfRepetitions' specifies the number of
%                   repetitions of the delta block (default =  1).
%
%               4.) 'InputChannel' specifies which input-channels of the
%                   uncertain plant are affected by delta. This should be
%                   specified as a cell:
%
%                                     [Rep1_Chnl_x,...,Rep1_Chnl_y]
%                     InputChannel ={ [    ...    ,...,    ...    ] }
%                                     [RepN_Chnl_x,...,RepN_Chnl_y]
%
%                    Here the order of the channels is not relevant.
%
%               4.) 'OutputChannel' specifies which output-channels of the
%                   uncertain plant are affected by delta. This should be
%                   specified as:
%
%                                     [Rep1_Chnl_x,...,Rep1_Chnl_y]
%                    OutputChannel ={ [    ...    ,...,    ...    ] }
%                                     [RepN_Chnl_x,...,RepN_Chnl_y]
%
%                    Here the order of the channels is not relevant.
%
%              In addition, various properties of the corresponding
%              IQC-multiplier can be specified: 
%
%               5.) 'BasisFunctionType' specifies the type of basis
%                   function (default = 1).
%
%               6.) 'Length' specifies the length of the basis function
%                   (default =  1).
%                                  
%               7.) 'PoleLocation' specifies the pole location of the basis
%                   function (default = -1).
%
%               8.) 'SampleTime' specifies the sampling time of the plant
%                   (default = 0, is continuous time)
%
%               9.) 'TerminalCost' defines the terminal cost constraint, Z,
%                   which is used in invariance analysis. Options are
%
%                    a.) 'off' (default)
%                    a.) 'on'
%
%              10.) 'PrimalDual' specifies whether the multiplier should be
%                   a primal/dual parametrization (default = 'Primal').
%                   
%                   Note: For a standard IQC-analysis, all multipliers are
%                   Primal ones.
% -------------------------------------------------------------------------

properties (SetAccess = public)
    % Uncertainty name
    Name string = 'delta';
    
    % Basis function type
    BasisFunctionType double {mustBeReal,mustBeFinite} = 1;

    % Length of the muliplier basis function
    Length double {mustBeReal,mustBeFinite} = 1;

    % Pole location of the multiplier basis function
    PoleLocation double {mustBeReal,mustBeFinite} = -1;
    
    % Sampling time of the multiplier basis function
    SampleTime double {mustBeReal,mustBeFinite} = 0;

    % Number of repetitions of the uncertain parameter
    NumberOfRepetitions double {mustBeReal,mustBeFinite} = 1;

    % Size of the uncertainty delta
    Dimensions double {mustBeReal,mustBeFinite} = 1;
   
    % Norm alpha of the uncertainty delta
    NormBounds double {mustBeReal,mustBeFinite} = 1;

    % Terminal cost certificate Z
    TerminalCost string {mustBeMember(TerminalCost,{'on','off'})} = 'off';
    
    % Primal or Dual multiplier
    PrimalDual string {mustBeMember(PrimalDual,{'Primal','Dual'})} = 'Primal';
    
    % InputChannel
    InputChannel = {};
    
    % OutputChannel
    OutputChannel = {};
end
methods
    function obj = ultid(varargin)
        if nargin == 0
            error('Error: To utilize "ultid" one should specify at least one input argument (see "help ultid" for further details)');
        elseif nargin == 1
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help ultid" for further details)');
            end
        elseif nargin > 1
            if nargin == 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help ultid" for further details)');
                end
                if isstruct(varargin{2})
                else
                    error('Error: The second input argument should be defined by a structure (see "help ultid" for further details)');
                end
                if isfield(varargin{2},'BasisFunctionType')
                    a = 1:1:10;
                    if isreal(varargin{2}.BasisFunctionType) && isscalar(varargin{2}.BasisFunctionType) && ismember(varargin{2}.BasisFunctionType,a)
                        obj.BasisFunctionType = varargin{2}.BasisFunctionType;
                    else
                        error('Error: The type of basis function should be defined as Natural numbers in N1');
                    end
                end
                if isfield(varargin{2},'Length')
                    a = 1:1:10;
                    if isreal(varargin{2}.Length) && isscalar(varargin{2}.Length) && ismember(varargin{2}.Length,a)
                        obj.Length = varargin{2}.Length;
                    else
                        error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultid" for further details).');
                    end
                end
                if isfield(varargin{2},'PoleLocation')
                    if isreal(varargin{2}.PoleLocation) && isscalar(varargin{2}.PoleLocation)
                        obj.PoleLocation = varargin{2}.PoleLocation;
                    else
                        error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultid" for further details).');
                    end
                end
                if isfield(varargin{2},'SampleTime')
                    if isreal(varargin{2}.SampleTime) && ismatrix(varargin{2}.SampleTime)
                        obj.SampleTime = varargin{2}.SampleTime;
                    else
                        error('Error: The Sampling time of a basis function should be defined as a nonnegative real rational number (see "help ultid" for further details).');
                    end
                end
                if isfield(varargin{2},'NumberOfRepetitions')
                    if isreal(varargin{2}.NumberOfRepetitions) && isscalar(varargin{2}.NumberOfRepetitions) && varargin{2}.NumberOfRepetitions > 0
                        obj.NumberOfRepetitions = varargin{2}.NumberOfRepetitions;
                    else
                        error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help ultid" for further details).');
                    end
                end
                if isfield(varargin{2},'Dimensions')
                    a = varargin{2}.Dimensions;
                    if isreal(a) && size(a,1) == 1 && size(a,2) == 2 && a(1,1) > 0 && a(1,2) > 0
                        obj.Dimensions = varargin{2}.Dimensions;
                    else
                        error('Error: The size of the uncertainty should be defined as [n_out,n_in], with n_out and n_in being positive natural numbers in N1 (see "help ultid" for further details).');
                    end
                end
                if isfield(varargin{2},'NormBounds')
                    if isreal(varargin{2}.NormBounds) && isscalar(varargin{2}.NormBounds) && varargin{2}.NormBounds > 0
                        obj.NormBounds = varargin{2}.NormBounds;
                    else
                        error('Error: The norm of the uncertainty should be defined as a positive real rational number (see "help ultid" for further details).');
                    end
                end
                if isfield(varargin{2},'TerminalCost')
                    obj.TerminalCost   = varargin{2}.TerminalCost;
                end
                if isfield(varargin{2},'PrimalDual')
                    obj.PrimalDual = varargin{2}.PrimalDual;
                end
                if isfield(varargin{2},'InputChannel')
                    obj.InputChannel = varargin{2}.InputChannel;
                end
                if isfield(varargin{2},'OutputChannel')
                    obj.OutputChannel = varargin{2}.OutputChannel;
                end
            elseif nargin > 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help ultid" for further details)');
                end
                n = length(varargin);
                if mod(n-1,2) ~= 0
                    error('Error: The input arguments should be defined in pairs')
                else
                    j = linspace(3,n,n/2);
                end
                for i = 1:length(j)
                    prop = varargin{j(i)-1};
                    switch prop
                        case 'BasisFunctionType'
                            a = 1:1:10;
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && ismember(varargin{j(i)},a)
                                obj.BasisFunctionType = varargin{j(i)};
                            else
                                error('Error: The type of basis function should be defined as natural numbers in N1 (see "help ultid" for further details)');
                            end
                        case 'Length'
                            a = 1:1:10;
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && ismember(varargin{j(i)},a)
                                obj.Length = varargin{j(i)};
                            else
                                error('Error: The length of a basis function should be defined as a natural number in N1 (see "help ultid" for further details).');
                            end
                        case 'PoleLocation'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)})
                                obj.PoleLocation = varargin{j(i)};
                            else
                                error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help ultid" for further details).');
                            end
                        case 'SampleTime'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.SampleTime = varargin{j(i)};
                            else
                                error('Error: The sampling time of a basis function should be defined as a nonnegative real rational number (see "help ultid" for further details).');
                            end    
                        case 'NumberOfRepetitions'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                                obj.NumberOfRepetitions = varargin{j(i)};
                            else
                                error('Error: The number of repetitions of the uncertainty should be defined as a vector with positive natural numbers in N1 (see "help ultid" for further details).');
                            end
                        case 'Dimensions'
                            a = varargin{j(i)};
                            if isreal(a) && size(a,1) == 1 && size(a,2) == 2 && a(1,1) > 0 && a(1,2) > 0
                                obj.Dimensions = varargin{j(i)};
                            else
                                error('Error: The size of the uncertainty should be defined as [n_out,n_in], with n_out and n_in being positive natural numbers in N1 (see "help ultid" for further details).');
                            end
                        case 'NormBounds'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                                obj.NormBounds = varargin{j(i)};
                            else
                                error('Error: The norm of the uncertainty should be defined as a positive real rational number (see "help ultid" for further details).');
                            end
                        case 'TerminalCost'
                            obj.TerminalCost = varargin{j(i)};
                        case 'PrimalDual'
                            obj.PrimalDual = varargin{j(i)};
                        case 'InputChannel'
                            obj.InputChannel = varargin{j(i)};
                        case 'OutputChannel'
                            obj.OutputChannel = varargin{j(i)};
                    end
                end    
            end
        end
    end
    function prob = iqcltid_p(Delta,prob)
    % ---------------------------------------------------------------------
    % This function generates two possible IQC-multipliers for the class
    % ultid:
    %
    % 1.) Dyn. D-scalings for diagonally repeated scalar uncertainties
    % 2.) Dyn. D-scalings for diagonally repeated full-block uncertainties
    %
    % These are precisely implemented as discussed in Section 5.1 of: 
    %
    % [1] J. Veenman, C.W. Scherer, H. Koroglu, Robust stability and
    %     performance analysis based on integral quadratic constraints,
    %     European Journal of Control, Vol. 31, pp. 1-32, 2016.
    % ---------------------------------------------------------------------
    % Dynamic D-scalings for diagonally repeated scalar uncertainties
    % ---------------------------------------------------------------------
    %
    % For this type of uncertainties this function implements the following
    % multiplier parametrization:
    %
    %     [x]^*[ D, 0 ][ bet*kron(Phi,Inr), 0             ]
    %     [x]  [ 0, -D][ 0                 , kron(Phi,Inr)]
    %
    % with:
    % - kron(Psi,I_nr)^*D*kron(Psi,I_nr) > 0
    % - nr is the number of repetitions of the uncertainty
    % - Phi\in\RHinf is the basis function 
    %
    % ---------------------------------------------------------------------
    % Dynamic D-scalings for diagonally repeated full-block uncertainties
    % ---------------------------------------------------------------------
    %
    % For this type of uncertainties this function implements the following
    % multiplier parametrization:
    %
    %     [x]^*[ kron(D,I_ni), 0         ][ bet*kron(Psi,I_ni), 0         ]
    %     [x]  [ 0        , -kron(D,I_no)][ 0             , kron(Psi,I_no)]
    %
    % with:
    % - Psi^*D*Psi >= 0
    % - Psi\in\RHinf is the basis function 
    %
    % ---------------------------------------------------------------------

    % Define outer-factor (basis function) of IQC-multiplier
    bft              = get(Delta,'BasisFunctionType');
    nr               = get(Delta,'NumberOfRepetitions');
    l                = get(Delta,'Length');
    pl               = get(Delta,'PoleLocation');
    Ts               = get(Delta,'SampleTime');
    dimu             = get(Delta,'Dimensions');
    bet              = get(Delta,'NormBounds');
    tc               = get(Delta,'TerminalCost');
    lnr              = length(nr);
    ll               = length(l);
    lpl              = length(pl);
    
    if lnr > 1 || ll > 1 || lpl > 1
        error('Error: The properties "Length", "PoleLocation", and "NumberOfRepetitions" should be defined as scalar inputs.');
    end
    if Ts == 0 && pl == 0
        disp('Warning: In continuous time it is not allowed to have poles on the imaginary axis. The "PoleLocation" was set to its nominal value (i.e. -1).');
        pl           = -1;
    elseif Ts > 0 && abs(pl) == 1
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl           = 0;
    elseif Ts == - 1 && abs(pl) == 1
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl           = 0;
    end
    
    % Define KYP certificate structure (for continuous or discrete time analysis)
    if Ts == 0
        kyp          = [0,1;1,0];
    else
        kyp          = [-1,0;0,1];
    end
    if dimu(1,1) == 1 && dimu(1,2) == 1
        % Define basis functions
        Phi          = fBasis(l,pl,nr,bft,Ts);
        Phi1         = Phi;
        Phi2         = Phi;
        prob         = fPsi(prob,Phi1,Phi2);

        % Define LMI variables
        D            = iqcvar(prob,[size(Phi.c,1),size(Phi.c,1)],'symmetric');
        X            = iqcvar(prob,[size(Phi.a,1),size(Phi.a,1)],'symmetric');
        V            = blkdiag(kron(kyp,X),D);
        P            = blkdiag(D,-D);
        sc           = blkdiag(bet*eye(D.Dim(1)),eye(D.Dim(1)));
        
        % Define LMI constraints
        A            = fss2m(Phi);
        prob         = iqclmi(prob,V,1,0,A);
        switch tc
            case 'on'
                Z    = blkdiag(X,-X);
                Tz   = blkdiag(bet*eye(X.Dim(1)),eye(X.Dim(1)));
                prob = fZ(prob,Z,Tz);
        end
    else
        % Define basis functions
        Phi          = fBasis(l,pl,1,bft,Ts);
        Phi1         = fsskron(Phi,dimu(2)*nr);
        Phi2         = fsskron(Phi,dimu(1)*nr);
        prob         = fPsi(prob,Phi1,Phi2);
        
        % Define LMI variables
        D            = iqcvar(prob,[size(Phi.c,1),size(Phi.c,1)],'symmetric');
        X            = iqcvar(prob,[size(Phi.a,1),size(Phi.a,1)],'symmetric');
        V            = blkdiag(kron(kyp,X),D);
        P11          = kron(D,eye(dimu(2)*nr));
        P22          = kron(-D,eye(dimu(1)*nr));
        P            = blkdiag(P11,P22);
        sc           = blkdiag(bet*eye(P11.Dim(1)),eye(P22.Dim(1)));
        
        % Define LMI constraints
        A            = fss2m(Phi);
        prob         = iqclmi(prob,V,1,0,A);
        switch tc
            case 'on'
                Z11  = kron(X,eye(dimu(2)*nr));
                Z22  = kron(-X,eye(dimu(1)*nr));
                Z    = blkdiag(Z11,Z22);
                Tz   = kron(blkdiag(bet*eye(dimu(2)*nr),eye(dimu(1)*nr)),eye(X.Dim(1)));
                prob = fZ(prob,Z,Tz);
        end
    end
    IO               = [dimu(2)*nr,dimu(1)*nr,0,0,size(Phi1.c,1),size(Phi2.c,1),size(Phi1.c,1),size(Phi2.c,1),size(Phi1.a,1),size(Phi2.a,1)];
    prob             = fP(prob,P);
    prob             = fsc(prob,sc);
    prob             = fIO(prob,IO);
    end
    function prob = iqcltid_d(Delta,prob)
        error('Error: Dual IQC-multipliers are not supported for this class of uncertainties.');
    end
    function display(obj)
        a1 = [];b1 = [];
        for i = 1:length(obj.OutputChannel)
            if i == 1
                a2 = ['[',num2str(obj.InputChannel{i}),']'];
                b2 = ['[',num2str(obj.OutputChannel{i}),']'];
            else
                a2 = [', [',num2str(obj.InputChannel{i}),']'];
                b2 = [', [',num2str(obj.OutputChannel{i}),']'];
            end
               a1 = [a1,a2];
               b1 = [b1,b2];
        end
        disp('--------------------------------------------------');
        disp(' LTI dynamic diagonally repeated uncertainty block');
        disp('--------------------------------------------------');
        disp(' Uncertainty details:');
        disp(['  - Name:                         ',num2str(obj.Name)]);
        disp(['  - Number of repetitions:        ',num2str(obj.NumberOfRepetitions)]);
        disp(['  - Input channels plant:         ',a1]);
        disp(['  - Output channels plant:        ',b1]);
        disp(['  - Dimension:                    ',num2str(obj.Dimensions(1,2)),' x ',num2str(obj.Dimensions(1,1))]);
        disp(['  - Norm uncertainty:             ',num2str(obj.NormBounds)]);
        disp(' ');
        disp(' Multiplier details:');
        disp(['  - Basis function type:          ',num2str(obj.BasisFunctionType)]);
        disp(['  - Basis function length:        ',num2str(obj.Length)]);
        disp(['  - Basis function pole location: ',num2str(obj.PoleLocation)]);
        disp(['  - Multiplier type:              ',num2str(obj.PrimalDual)]);
    end
end 
end