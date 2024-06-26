classdef usg < handle & matlab.mixin.SetGetExactNames

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
% Date:        02-12-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition for norm-bounded uncertainties Delta:L2->L2
%              with: 
%              - Dimension n_out x n_in
%              - ||Delta||_L2 < alpha
%
% Syntax:      Delta = usg('name')
%              Delta = usg('name',varargin)
%
% Usage:       "Delta = usg('name')" defines a SISO (possibly nonlinear)
%              norm-bounded uncertainty with L2-norm equal to 1.
%
%              For "Delta = usg('name',varargin)", the varargin inputs come
%              in pairs and can be defined as: 
%
%              delta = usg('name','prop1','value1','prop2','value2',...)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              usgOpt.prop1 = value1
%              usgOpt.prop2 = value2
%                       ...              
%              usgOpt.propN = valueN
% 
%              Delta = usg('name',usgOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%              set(Delta,'propX','valueX') and get(Delta,'propX').
%
%              The properties that can be specified are:
%            
%               1.) 'Dimensions' specifies the size of the uncertainty
%                   delta as [n_out, n_in] (default = [1,1]).
%
%               2.) 'NormBounds' specifies the L2-norm of the uncertainty
%                   by means of the value alpha (default = 1).
%
%               3.) 'InputChannel' specifies which input-channels of the
%                   uncertain plant are affected by Delta. The channels
%                   should be specified as a cell:
%
%                           InputChannel = {[Chnl_x,...,Chnl_y]}
%
%                    Here the order of the channels is not relevant, while
%                    the row-length equals the number of outputs of Delta.
%
%               4.) 'OutputChannel' specifies which input-channels of the
%                   uncertain plant are affected by Delta. The channels
%                   should be specified as a cell:
%
%                           OutputChannel = {[Chnl_x,...,Chnl_y]}
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
%                   Note: The length of the basis function can only be
%                   larger than 1 if the integrant
%
%                       [I-Delta(q(i\omega))^*Delta(q(i\omega))]\geq0 
%
%                   for all omega\in R\cup{\infty}.
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
    Name string = 'Delta';
    
    % Basis function type
    BasisFunctionType double {mustBeReal,mustBeFinite} = 1;

    % Length of the muliplier basis function
    Length double {mustBeReal,mustBeFinite} = 1;

    % Pole location of the multiplier basis function
    PoleLocation double {mustBeReal,mustBeFinite} = -1;
    
    % Sampling time of the multiplier basis function
    SampleTime double {mustBeReal,mustBeFinite} = 0;

    % Size of the uncertainty delta
    Dimensions double {mustBeReal,mustBeFinite} = [1,1];
   
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
    function obj = usg(varargin)
        if nargin == 0
            error('Error: To utilize "usg" one should specify at least one input argument (see "help usg" for further details)');
        elseif nargin == 1
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help usg" for further details)');
            end
        elseif nargin > 1
            if nargin == 2
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('Error: The name of the uncertainty should be defined by a string (see "help usg" for further details)');
                end
                if isstruct(varargin{2})
                else
                    error('Error: The second input argument should be defined by a structure (see "help usg" for further details)');
                end
                if isfield(varargin{2},'BasisFunctionType')
                    a = 1:1:5;
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
                        if varargin{2}.Length > 1
                            disp('Warning: "Length" is only allowed to be larger than 0 if the integrant of the IQC is positive for all omega in [0,\infty].');
                        end
                    else
                        error('Error: The length of a basis function should be defined as a natural number in N1 (see "help usg" for further details).');
                    end
                end
                if isfield(varargin{2},'PoleLocation')
                    if isreal(varargin{2}.PoleLocation) && isscalar(varargin{2}.PoleLocation)
                        obj.PoleLocation = varargin{2}.PoleLocation;
                    else
                        error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help usg" for further details).');
                    end
                end
                if isfield(varargin{2},'SampleTime')
                    if isreal(varargin{2}.SampleTime) && ismatrix(varargin{2}.SampleTime)
                        obj.SampleTime = varargin{2}.SampleTime;
                    else
                        error('Error: The Sampling time of a basis function should be defined as a nonnegative real rational number (see "help usg" for further details).');
                    end
                end
                if isfield(varargin{2},'Dimensions')
                    a = varargin{2}.Dimensions;
                    if isreal(a) && size(a,1) == 1 && size(a,2) == 2 && a(1,1) > 0 && a(1,2) > 0
                        obj.Dimensions = varargin{2}.Dimensions;
                    else
                        error('Error: The size of the uncertainty should be defined as [n_out,n_in], with n_out and n_in being positive natural numbers in N1 (see "help usg" for further details).');
                    end
                end
                if isfield(varargin{2},'NormBounds')
                    if isreal(varargin{2}.NormBounds) && isscalar(varargin{2}.NormBounds) && varargin{2}.NormBounds > 0
                        obj.NormBounds = varargin{2}.NormBounds;
                    else
                        error('Error: The norm of the uncertainty should be defined as a positive real rational number (see "help usg" for further details).');
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
                    error('Error: The name of the uncertainty should be defined by a string (see "help usg" for further details)');
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
                            a = 1:1:5;
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && ismember(varargin{j(i)},a)
                                obj.BasisFunctionType = varargin{j(i)};
                            else
                                error('Error: The type of basis function should be defined as natural numbers in N1 (see "help usg" for further details)');
                            end
                        case 'Length'
                            a = 1:1:10;
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && ismember(varargin{j(i)},a)
                                obj.Length = varargin{j(i)};
                                if varargin{j(i)} > 1
                                    disp('Warning: "Length" is only allowed to be larger than 0 if the integrant of the IQC is positive for all omega in [0,\infty].');
                                end
                            else
                                error('Error: The length of a basis function should be defined as a natural number in N1 (see "help usg" for further details).');
                            end
                        case 'PoleLocation'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)})
                                obj.PoleLocation = varargin{j(i)};
                            else
                                error('Error: The pole location of a basis function should be defined as a real rational number. For continuous and distrete time systems, it is not allowed to have poles on the imaginary axis and unit circle respectively (see "help usg" for further details).');
                            end
                        case 'SampleTime'
                            if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                                obj.SampleTime = varargin{j(i)};
                            else
                                error('Error: The sampling time of a basis function should be defined as a nonnegative real rational number (see "help usg" for further details).');
                            end    
                        case 'Dimensions'
                            a = varargin{j(i)};
                            if isreal(a) && size(a,1) == 1 && size(a,2) == 2 && a(1,1) > 0 && a(1,2) > 0
                                obj.Dimensions = varargin{j(i)};
                            else
                                error('Error: The size of the uncertainty should be defined as [n_out,n_in], with n_out and n_in being positive natural numbers in N1 (see "help usg" for further details).');
                            end
                        case 'NormBounds'
                            if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                                obj.NormBounds = varargin{j(i)};
                            else
                                error('Error: The norm of the uncertainty should be defined as a positive real rational number (see "help usg" for further details).');
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
    function prob = iqcsg_p(Delta,prob)
    % ---------------------------------------------------------------------
    % This function generates an IQC-multipliers for the class of
    % full-block norm-bounded nonlinearities. See [1], Section 5.7 for
    % further details:
    %
    % [1] J. Veenman, C.W. Scherer, H. Koroglu, Robust stability and
    %     performance analysis based on integral quadratic constraints,
    %     European Journal of Control, Vol. 31, pp. 1-32, 2016.
    %
    % ---------------------------------------------------------------------

    % Define outer-factor (basis function) of IQC-multiplier
    bft           = get(Delta,'BasisFunctionType');
    l             = get(Delta,'Length');
    pl            = get(Delta,'PoleLocation');
    Ts            = get(Delta,'SampleTime');
    dimu          = get(Delta,'Dimensions');
    bet           = get(Delta,'NormBounds');
    tc            = get(Delta,'TerminalCost');
    ll            = length(l);
    lpl           = length(pl);
    
    if l > 1
        disp('Warning: The property "Length" can only be set to values larger than 1 if the integrant I-Delta(q(i\omega))^*Delta(q(i\omega)) is greater or equal to zero for all omega\in R\cup{\infty}.');
    end
    if ll > 1 || lpl > 1
        error('Error: The properties "Length", "PoleLocation", and "NumberOfRepetitions" should be defined as scalar inputs.');
    end
    if Ts == 0 && pl == 0
        disp('Warning: In continuous time it is not allowed to have poles on the imaginary axis. The "PoleLocation" was set to its nominal value (i.e. -1).');
        pl        = -1;
    elseif Ts > 0 && abs(pl) == 1
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl        = 0;
    elseif Ts == -1 && abs(pl) == 1
        disp('Warning: In discrete time it is not allowed to have poles on the unit circle. The "PoleLocation" was set to its nominal value (i.e. 0).');
        pl        = 0;
    end

    % Define KYP certificate structure (for continuous or discrete time analysis)
    if Ts == 0
        kyp       = [0,1;1,0];
    else
        kyp       = [-1,0;0,1];
    end
    
    % Define basis functions
    Phi           = fBasis(l,pl,1,bft,Ts);
    Phi1          = fsskron(Phi,dimu(2));
    Phi2          = fsskron(Phi,dimu(1));
    prob          = fPsi(prob,Phi1,Phi2);
    
    % Define LMI variables
    D             = iqcvar(prob,[size(Phi.c,1),size(Phi.c,1)],'symmetric');
    X             = iqcvar(prob,[size(Phi.a,1),size(Phi.a,1)],'symmetric');
    V             = blkdiag(kron(kyp,X),D);
    P11           = kron(D,eye(dimu(2)));
    P22           = kron(-D,eye(dimu(1)));
    P             = blkdiag(P11,P22);
    sc            = blkdiag(bet*eye(P11.Dim(1)),eye(P22.Dim(1)));
    prob          = fP(prob,P);
    prob          = fsc(prob,sc);
    
    % Define LMI constraints
    A             = fss2m(Phi);
    prob          = iqclmi(prob,V,1,0,A);
    
    switch tc
        case 'on'
            xPhi  = size(Phi1.a,1) + size(Phi2.a,1);
            Z     = iqcvar(prob,[xPhi,xPhi],'symmetric');
            prob  = fZ(prob,Z,eye(Z.Dim(1)));
            nPhi1 = size(Phi1,2);
            nPhi2 = size(Phi2,2);
            B     = fss2m(blkdiag(Phi1,Phi2)*[zeros(nPhi1,nPhi2);eye(nPhi2)]);
            W     = blkdiag(kron(kyp,Z),P);
            prob  = iqclmi(prob,W,-1,0,B);
    end
    
    % Define IO data of the multipliers
    IO            = [dimu(2),dimu(1),0,0,size(Phi1.c,1),size(Phi2.c,1),size(Phi1.c,1),size(Phi2.c,1),size(Phi1.a,1),size(Phi2.a,1)];
    prob          = fIO(prob,IO);
    end
    function prob = iqcsg_d(Delta,prob)
        error('Error: Dual IQC-multipliers are not supported for this class of uncertainties.');
    end
    function display(obj)
        a1    = [];b1 = [];
        for i = 1:length(obj.OutputChannel)
            if i == 1
                a2 = ['[',num2str(obj.InputChannel{i}),']'];
                b2 = ['[',num2str(obj.OutputChannel{i}),']'];
            else
                a2 = [', [',num2str(obj.InputChannel{i}),']'];
                b2 = [', [',num2str(obj.OutputChannel{i}),']'];
            end
               a1  = [a1,a2];
               b1  = [b1,b2];
        end
        lout       = length(obj.OutputChannel{1});
        lin        = length(obj.InputChannel{1});
        dimDelta   = [lin,lout];
        disp('-----------------------------------------');
        disp(' Nonlinear norm-bounded uncertainty block');
        disp('-----------------------------------------');
        disp(' Uncertainty details:');
        disp(['  - Name:                         ',num2str(obj.Name)]);
        disp(['  - Input channels plant:         ',a1]);
        disp(['  - Output channels plant:        ',b1]);
        disp(['  - Dimension:                    ',num2str(dimDelta(1,1)),' x ',num2str(dimDelta(1,2))]);
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