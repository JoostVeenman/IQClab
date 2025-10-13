classdef iqcdelta < matlab.mixin.SetGetExactNames

% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        29-11-2019
% 
% -------------------------------------------------------------------------
%
% Description: Class definition of uncertainty and performance blocks
%
% Syntax:      delta = iqcdelta('name')
%              delta = iqcdelta('name',varargin)
%
% Usage:       "delta = iqcdelta('name')" defines an empty uncertainty. For
%              "delta = iqcdelta('name',varargin)", the varargin inputs 
%              come in pairs and can be defined as: 
%
%              delta = iqcdelta('name','prop1','value1','prop2','value2',...)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              iqcdeltaOpt.prop1 = value1
%              iqcdeltaOpt.prop2 = value2
%                         ...              
%              iqcdeltaOpt.propN = valueN
% 
%              delta = iqcdelta('name',iqcdeltaOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%              set(delta,'propX','valueX') and get(delta,'propX').
%
%              It is assumed that the uncertain plant LFT(Delta,M) is known
%              where M is the uncertain LTI plant and Delta the overall
%              uncertainty block.
%
%              The properties that can be specified are:
%
%          ----------------------------------------------------------------
%               A.) Properties related to IO-channels of plant M
%          ----------------------------------------------------------------
%
%               1.) 'ChannelClass' specifies whether the considered
%                   channels are uncertainty channels 'U',performance
%                   channels 'P', or control channels 'C' (default: 'U').
%               2.) 'InputChannel' and 'OutputChannel' specify the input
%                   and output channels of the plant M:
%
%                   a.) For a single repeated SISO uncertainty block:
%
%                       InputChannel  = [chnlx,...,chnly]
%                       OutputChannel = [chnlv,...,chnlw]
%
%                       Note: the order of the channels is not important
%
%                   b.) For a single repeated MIMO uncertainty block:
%
%                                       [chnl_delta_x1, ... ,chnl_delta_y1]
%                       InputChannel  = [      ...    , ... ,     ...     ]
%                                       [chnl_delta_xN, ... ,chnl_delta_yN]
%
%                                       [chnl_delta_v1, ... ,chnl_delta_v1]
%                       OutputChannel = [      ...    , ... ,     ...     ]
%                                       [chnl_delta_wN, ... ,chnl_delta_wN]
%
%                       Note: the order of the channels is not important
%
%                   c.) For multiple repeated SISO or MIMO uncertainty
%                       blocks, one should proceed as before with the
%                       difference that each row should be specified as the
%                       element of a cell:
%
%                       InputChannel  = {row1 , ... , rowN}
%                       OutputChannel = {row1 , ... , rowN}
%
%          ----------------------------------------------------------------
%               B.) Properties related to the nature of uncertainty
%          ----------------------------------------------------------------
%
%               1.) 'LinNonlin' spefifies whether the uncertainty is linear
%                   ('L') or nonlinear ('NL') - (default: 'L').
%               2.) 'TimeInvTimeVar' specifies whether the uncertainty is
%                   time invariant ('TI') or time varying ('TV') -
%                   (default: 'TI'). 
%               3.) 'StaticDynamic' specified whether the uncertainty is
%                   static ('S') or dynamic ('D') - (default: 'S').
%               4.) 'Structure' specifies whether the uncertainty is
%                   diagonal ('D') or full block ('FB') - (default: 'D').
%
%          ----------------------------------------------------------------
%               C.) Properties related to the characteristics of uncertainty
%          ----------------------------------------------------------------
%
%               1.) 'Bounds' (if any) specifies the bounds on an uncertain
%                    parameter:
%
%                    a.) For a single parameter:
%
%                        Bounds = [-b,b]
%
%                    b.) For multiple parameters:
%
%                        Bounds = {[-b1,b1], ... , [-bN,bN]}
%
%               2.) 'RateBounds' (if any) specifies the rate bounds on an
%                    uncertain parameter:
%
%                    a.) For a single parameter:
%
%                        RateBounds = [-rb,rb]
%
%                    b.) For multiple parameters:
%
%                        RateBounds = {[-rb1,rb1], ... , [-rbN,rbN]}
%
%               3.) 'NormBounds' (if any) specifies the norm bounds on an
%                    uncertain parameter:
%
%                    a.) NormBounds = nb
%
%               4.) 'SectorBounds' (if any) specifies the sector bounds on
%                    an uncertainty/nonlinearity:
%
%                    a.) For a single uncertainty/nonlinearity:
%
%                        SectorBounds = [a,b], a<=0<=b
%
%                    b.) For multiple uncertainty/nonlinearity:
%
%                        SectorBounds = {[a1,b1], ... , [aN,bN]}, ai<=0<=bi
%
%               5.) 'SlopeBounds' (if any) specifies the slope bounds on
%                    an uncertainty/nonlinearity:
%
%                    a.) For a single uncertainty/nonlinearity:
%
%                        SlopeBounds = [a,b], a<=0<=b
%
%                    b.) For multiple uncertainty/nonlinearity:
%
%                        SlopeBounds = {[a1,b1], ... , [aN,bN]}, ai<=0<=bi
%
%               6.) 'Polytope' specifies the generator points of an
%                    uncertainty block to be specified as:
%
%                                   [delta1^1, ... ,deltaN^1]
%                    a.) Polytope = [   ...  , ... ,   ...  ]
%                                   [delta1^M, ... ,deltaN^M]
%
%               7.) 'UncertaintyMap' specifies the map of an uncertainty by
%                    means of the matrices T1,...,TN as a cell array T =
%                    {T1,...,TN}. This defines the uncertainty mapping: 
%
%                    a.) Delta(delta) = sum_{i=1}^N delta_i*Ti = ...
%                                       = delta1*T1 + ... + deltaN*TN.
%
%               8.) 'DelayTime' specifies the maximum delay time of a delay
%                    uncertainty:
%
%                    a.) DelayTime = value
%
%               9.) 'DelayType' specifies the type of delay operator (see
%                    udel for further details):
%
%                    a.) DelayType = value
%
%              10.) 'Passive' specifies whether an uncertainty block is
%                    passive 
%
%                    a.) 'Passive'
%
%              11.) 'Odd' specifies whether a nonlinear block is odd (see
%                    usbsr for further details). 
%
%                    Odd = 'yes' (i.e.  f(-x) = -f(x)) (else Odd = 'no')
%
%              12.) 'PerfMetric' specifies the particular performance
%                    metric for the performance channels.
%
%                    a.) PerfMetric = 'L2', 'H2', 'genH2', 'Passive',
%                                     'e2x', 'e2p', 'x2z', 'e2z', 'x2p'
%
%          ----------------------------------------------------------------
%               D.) Combining uncertainty blocks
%          ----------------------------------------------------------------
%                   It is possible to diagonally combine uncertainty blocks
%                   by means of the function
%
%                   Delta = combinedeltas('name',block1,...,blockN)
%
%                   The function combines all properties assigned to the
%                   individual blocks to create a new and larger one which
%                   carries all properties of the individual blocks.
%
% -------------------------------------------------------------------------

% Properties related to the channel selection
properties  (SetAccess = private)
    Name           string                                                                                     = 'delta'; % Name uncertainty
    ChannelClass   string {mustBeMember(ChannelClass,{'U','P','C'})}                                          = 'U';     % Uncertainty (U), Performance (P), or Control (C) channel
    InputChannel   cell                                                                                       = {};      % Input channel selection plant M of LFT(Delta,M)
    OutputChannel  cell                                                                                       = {};      % Output channel selection plant M of LFT(Delta,M)
end

% Properties related to the nature of the uncertainty/performance channels
properties (SetAccess = private)   
    LinNonlin      cell                                                                                       = {'L'};   % Linear (L) or NonLinear (NL)
    TimeInvTimeVar cell                                                                                       = {'TI'};  % Time-Invariant (TI) or Time-Varying (TV)
    StaticDynamic  cell                                                                                       = {'S'};   % Static (S) or Dynamnic (D)
    Structure      cell                                                                                       = {'D'};   % Diagonal (D) or FullBlock (FB)
end
% Properties related to the characteristics of uncertainties/performances
properties (SetAccess = public)    
    Bounds         cell                                                                                       = {[]};    % Bounds on the uncertainty (interval should contain 0)
    RateBounds     cell                                                                                       = {[]};    % Rate bounds on the uncertainty  (interval should contain 0)
    NormBounds     cell                                                                                       = {[]};    % Norm bounds on the uncertainty
    SectorBounds   cell                                                                                       = {[]};    % Sector bounds of the uncertainty  (interval should contain 0)
    SlopeBounds    cell                                                                                       = {[]};    % Slope restriction of the uncertainty
    Polytope       double                                                                                     = [];      % Polytope defining uncertainty region (interval should contain 0)
    UncertaintyMap cell                                                                                       = {};      % Uncertainty Mapping
    DelayTime      cell                                                                                       = {[]};    % Delay time
    DelayType      cell                                                                                       = {[]};    % Delay type
    Passive        string {mustBeMember(Passive,{'Passive',''})}                                              = {''};    % Passiveness of the uncertainty
    Odd            string {mustBeMember(Odd,{'yes','no',''})}                                                 = {''};    % Specify if a nonlinearity is an odd or even function
    PerfMetric     string {mustBeMember(PerfMetric,{'L2','H2','genH2','Passive','e2x','e2p','x2z','e2z','x2p',''})} = {'L2'};  % Performance metric
end
methods
    function obj = iqcdelta(varargin)
        if nargin == 0
            error('Error: One should at least specify the name of the uncertainty');
        elseif nargin == 1
            if ischar(varargin{1})
                obj.Name = varargin{1};
            elseif isobject(varargin{1})
                obj = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help iqcdelta" for further details)');
            end
        elseif nargin == 2
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help iqcdelta" for further details)');
            end
            if isstruct(varargin{2})
            else
                error('Error: The second input argument should be defined by a structure (see "help iqcdelta" for further details)');
            end
            % Set channel properties
            if isfield(varargin{2},'ChannelClass')
                obj.ChannelClass          = varargin{2}.ChannelClass;
                if strcmp(varargin{2}.ChannelClass,'P')
                    obj.LinNonlin{1}      = '';
                    obj.TimeInvTimeVar{1} = '';
                    obj.StaticDynamic{1}  = '';
                    obj.Structure{1}      = '';
                elseif strcmp(varargin{2}.ChannelClass,'U')
                    obj.PerfMetric        = '';
                elseif strcmp(varargin{2}.ChannelClass,'C')
                    obj.LinNonlin{1}      = '';
                    obj.TimeInvTimeVar{1} = '';
                    obj.StaticDynamic{1}  = '';
                    obj.Structure{1}      = '';
                    obj.PerfMetric        = '';
                end
            end
            if isfield(varargin{2},'InputChannel')
                InCh   = varargin{2}.InputChannel;
                [s1,~] = size(InCh);
                if s1 == 1
                    if iscell(InCh)
                        for ijk = 1:length(InCh)
                            obj.InputChannel{ijk} = InCh{ijk};
                        end
                    else 
                        obj.InputChannel{1} = InCh;
                    end
                else
                    if iscell(InCh)
                        for ijk = 1:length(InCh)
                            obj.InputChannel{ijk} = InCh{ijk};
                        end
                    else 
                        for ijk = 1:s1
                            obj.InputChannel{ijk} = InCh(ijk,:);
                        end
                    end
                end
                if strcmp(obj.ChannelClass,'U')
                    obj.PerfMetric = '';
                end
            end
            if isfield(varargin{2},'OutputChannel')
                OutCh  = varargin{2}.OutputChannel;
                [s1,~] = size(OutCh);
                if s1 == 1
                    if iscell(OutCh)
                        for ijk = 1:length(OutCh)
                            obj.OutputChannel{ijk} = OutCh{ijk};
                        end
                    else 
                        obj.OutputChannel{1} = OutCh;
                    end
                else
                    if iscell(OutCh)
                        for ijk = 1:length(OutCh)
                            obj.OutputChannel{ijk} = OutCh{ijk};
                        end
                    else 
                        for ijk = 1:s1
                            obj.OutputChannel{ijk} = OutCh(ijk,:);
                        end
                    end
                end
                if strcmp(obj.ChannelClass,'U')
                    obj.PerfMetric = '';
                end
            end

            % Set uncertainty nature properties
            if isfield(varargin{2},'LinNonlin')
                s1                        = {'L','NL'};
                if ismember(varargin{2}.LinNonlin,s1)
                    obj.LinNonlin{1}      = varargin{2}.LinNonlin;
                else
                    error('Error: The property "LinNonlin", can only assume the characters "L" or "NL".');
                end
            end
            if isfield(varargin{2},'TimeInvTimeVar')
                s1                        = {'TI','TV'};
                if ismember(varargin{2}.TimeInvTimeVar,s1)
                    obj.TimeInvTimeVar{1} = varargin{2}.TimeInvTimeVar;
                else
                    error('Error: The property "TimeInvTimeVar", can only assume the characters "TI" or "TV".');
                end
            end
            if isfield(varargin{2},'StaticDynamic')
                s1                        = {'S','D'};
                if ismember(varargin{2}.StaticDynamic,s1)
                    obj.StaticDynamic{1}  = varargin{2}.StaticDynamic;
                else
                    error('Error: The property "StaticDynamic", can only assume the characters "S" or "D".');
                end
            end
            if isfield(varargin{2},'Structure')
                s1                        = {'D','FB'};
                if ismember(varargin{2}.Structure,s1)
                    obj.Structure{1}      = varargin{2}.Structure;
                else
                    error('Error: The property "Structure", can only assume the characters "D" or "FB".');
                end
            end
            
            % Set uncertainty characteristics
            if isfield(varargin{2},'Bounds')
                if iscell(varargin{2}.Bounds)
                    obj.Bounds            = varargin{2}.Bounds;
                else
                    obj.Bounds{1}         = varargin{2}.Bounds;
                end
            end
            if isfield(varargin{2},'RateBounds')
                if iscell(varargin{2}.RateBounds)
                    obj.RateBounds        = varargin{2}.RateBounds;
                else
                    obj.RateBounds{1}     = varargin{2}.RateBounds;
                end
            end
            if isfield(varargin{2},'NormBounds')
                if iscell(varargin{2}.NormBounds)
                    obj.NormBounds        = varargin{2}.NormBounds;
                else
                    obj.NormBounds{1}     = varargin{2}.NormBounds;
                end
            end
            if isfield(varargin{2},'SectorBounds')
                if iscell(varargin{2}.SectorBounds)
                    obj.SectorBounds      = varargin{2}.SectorBounds;
                else
                    obj.SectorBounds{1}   = varargin{2}.SectorBounds;
                end
            end
            if isfield(varargin{2},'SlopeBounds')
                if iscell(varargin{2}.SlopeBounds)
                    obj.SlopeBounds       = varargin{2}.SlopeBounds;
                else
                    obj.SlopeBounds{1}    = varargin{2}.SlopeBounds;
                end
            end
            if isfield(varargin{2},'Polytope')
                obj.Polytope              = varargin{2}.Polytope;
            end
            if isfield(varargin{2},'UncertaintyMap')
                obj.UncertaintyMap        = varargin{2}.UncertaintyMap;
            end
            if isfield(varargin{2},'DelayTime')
                if iscell(varargin{2}.DelayTime)
                    obj.DelayTime         = varargin{2}.DelayTime;
                else
                    obj.DelayTime{1}      = varargin{2}.DelayTime;
                end
            end
            if isfield(varargin{2},'DelayType')
                if iscell(varargin{2}.DelayType)
                    obj.DelayType         = varargin{2}.DelayType;
                else
                    obj.DelayType{1}      = varargin{2}.DelayType;
                end
            end
            if isfield(varargin{2},'Passive')
                s1                        = {'Passive',''};
                if ismember(varargin{2}.Structure,s1)
                    obj.Passive{1}        = varargin{2}.Passive;
                else
                    error('Error: The property "Passive", can only assume the character "Passive".');
                end
            end
            if isfield(varargin{2},'Odd')
                s1                        = {'yes','no',''};
                if ismember(varargin{2}.Structure,s1)
                    obj.Odd{1}            = varargin{2}.Odd;
                else
                    error('Error: The property "Odd", can only assume the character "yes" or "no".');
                end
            end
            if isfield(varargin{2},'PerfMetric')
                obj.PerfMetric            = varargin{2}.PerfMetric;
            end
        elseif nargin > 2 % specify properties via input-pairs
            if ischar(varargin{1})
                obj.Name = varargin{1};
            else
                error('Error: The name of the uncertainty should be defined by a string (see "help iqcdelta" for further details)');
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
                    % Set channel properties
                    case 'ChannelClass'
                        obj.ChannelClass          = varargin{j(i)};
                        if strcmp(varargin{j(i)},'P')
                            obj.LinNonlin{1}      = '';
                            obj.TimeInvTimeVar{1} = '';
                            obj.StaticDynamic{1}  = '';
                            obj.Structure{1}      = '';
                        elseif strcmp(varargin{j(i)},'U')
                            obj.PerfMetric        = '';
                        elseif strcmp(varargin{j(i)},'C')
                            obj.LinNonlin{1}      = '';
                            obj.TimeInvTimeVar{1} = '';
                            obj.StaticDynamic{1}  = '';
                            obj.Structure{1}      = '';
                            obj.PerfMetric        = '';
                        end
                    case 'InputChannel'
                        InCh   = varargin{j(i)};
                        [s1,~] = size(InCh);
                        if s1 == 1
                            if iscell(InCh)
                                for ijk = 1:length(InCh)
                                    obj.InputChannel{ijk} = InCh{ijk};
                                end
                            else 
                                obj.InputChannel{1}       = InCh;
                            end
                        else
                            if iscell(InCh)
                                for ijk = 1:length(InCh)
                                    obj.InputChannel{ijk} = InCh{ijk};
                                end
                            else 
                                for ijk = 1:s1
                                    obj.InputChannel{ijk} = InCh(ijk,:);
                                end
                            end
                        end
                        if strcmp(obj.ChannelClass,'U')
                            obj.PerfMetric = '';
                        end
                    case 'OutputChannel'   
                        OutCh  = varargin{j(i)};
                        [s1,~] = size(OutCh);
                        if s1 == 1
                            if iscell(OutCh)
                                for ijk = 1:length(OutCh)
                                    obj.OutputChannel{ijk} = OutCh{ijk};
                                end
                            else 
                                obj.OutputChannel{1}       = OutCh;
                            end
                        else
                            if iscell(OutCh)
                                for ijk = 1:length(OutCh)
                                    obj.OutputChannel{ijk} = OutCh{ijk};
                                end
                            else 
                                for ijk = 1:s1
                                    obj.OutputChannel{ijk} = OutCh(ijk,:);
                                end
                            end
                        end
                        if strcmp(obj.ChannelClass,'U')
                            obj.PerfMetric        = '';
                        end

                    % Set uncertainty nature properties
                    case 'LinNonlin'
                        s1                        = {'L','NL'};
                        if ismember(varargin{j(i)},s1)
                            obj.LinNonlin{1}      = varargin{j(i)};
                        else
                            error('Error: The property "LinNonlin", can only assume the characters "L" or "NL".');
                        end
                    case 'TimeInvTimeVar'
                        s1                        = {'TI','TV'};
                        if ismember(varargin{j(i)},s1)
                            obj.TimeInvTimeVar{1} = varargin{j(i)};
                        else
                            error('Error: The property "TimeInvTimeVar", can only assume the characters "TI" or "TV".');
                        end
                    case 'StaticDynamic'
                        s1                        = {'S','D'};
                        if ismember(varargin{j(i)},s1)
                            obj.StaticDynamic{1}  = varargin{j(i)};
                        else
                            error('Error: The property "StaticDynamic", can only assume the characters "S" or "D".');
                        end
                    case 'Structure'
                        s1                        = {'D','FB'};
                        if ismember(varargin{j(i)},s1)
                            obj.Structure{1}      = varargin{j(i)};
                        else
                            error('Error: The property "Structure", can only assume the characters "D" or "FB".');
                        end

                    % Set uncertainty characteristics
                    case 'Bounds'
                        if iscell(varargin{j(i)})
                            obj.Bounds            = varargin{j(i)};
                        else
                            obj.Bounds{1}         = varargin{j(i)};
                        end
                    case 'RateBounds'
                        if iscell(varargin{j(i)})
                            obj.RateBounds        = varargin{j(i)};
                        else
                            obj.RateBounds{1}     = varargin{j(i)};
                        end
                    case 'NormBounds'
                        if iscell(varargin{j(i)})
                            obj.NormBounds        = varargin{j(i)};
                        else
                            obj.NormBounds{1}     = varargin{j(i)};
                        end
                    case 'SectorBounds'
                        if iscell(varargin{j(i)})
                            obj.SectorBounds      = varargin{j(i)};
                        else
                            obj.SectorBounds{1}   = varargin{j(i)};
                        end
                    case 'SlopeBounds'
                        if iscell(varargin{j(i)})
                            obj.SlopeBounds       = varargin{j(i)};
                        else
                            obj.SlopeBounds{1}    = varargin{j(i)};
                        end
                    case 'Polytope'
                        obj.Polytope              = varargin{j(i)};
                    case 'UncertaintyMap'
                        obj.UncertaintyMap        = varargin{j(i)};
                    case 'DelayTime'
                        if iscell(varargin{j(i)})
                            obj.DelayTime         = varargin{j(i)};
                        else
                            obj.DelayTime{1}      = varargin{j(i)};
                        end
                    case 'DelayType'
                        if iscell(varargin{j(i)})
                            obj.DelayType         = varargin{j(i)};
                        else
                            obj.DelayType{1}      = varargin{j(i)};
                        end
                    case 'Passive'
                        s1                        = {'Passive',''};
                        if ismember(varargin{j(i)},s1)
                            obj.Passive{1}        = varargin{j(i)};
                        else
                            error('Error: The property "Passive", can only assume the character "Passive".');
                        end
                    case 'Odd'
                        s1                        = {'yes','no',''};
                        if ismember(varargin{j(i)},s1)
                            obj.Odd{1}        = varargin{j(i)};
                        else
                            error('Error: The property "Odd", can only assume the character "yes" or "no".');
                        end    
                    case 'PerfMetric'
                        obj.PerfMetric            = varargin{j(i)};
                end
            end    
        end
    end
    function obj = blkdiag(Name,varargin)
        obj = iqcdelta(Name);
        n   = length(varargin);
        if n < 2
           error('Error: It is necessary to provide at least two uncertainties to construct a combined one.');
        end
        for i = 1:n
            if  i == 1
                obj.ChannelClass       = varargin{i}.ChannelClass;
                obj.PerfMetric         = varargin{i}.PerfMetric;
                if strcmp(varargin{i}.ChannelClass,'P') || strcmp(varargin{i}.ChannelClass,'C')
                    obj.LinNonlin      = varargin{i}.LinNonlin;
                    obj.TimeInvTimeVar = varargin{i}.TimeInvTimeVar;
                    obj.StaticDynamic  = varargin{i}.StaticDynamic;
                    obj.Structure      = varargin{i}.Structure;
                end
            elseif i > 1 
                m = strcmp(varargin{i}.ChannelClass,varargin{i-1}.ChannelClass);
                k = strcmp(varargin{i}.PerfMetric,varargin{i-1}.PerfMetric);
                if m == 0
                    error('Error: One can not combine different channel classes.');
                end
                if k == 0
                    error('Error: One can not consider multi-objective performance criteria.');
                end
            end
            if isempty(varargin{i}.Polytope)
                obj.Polytope           = [];
                else
                    if isempty(obj.Polytope)
                        obj.Polytope   = varargin{i}.Polytope;
                    else
                        [sa1,~]        = size(obj.Polytope);
                        [sb1,~]        = size(varargin{i}.Polytope);
                        d              = kron(ones(sb1,1),obj.Polytope);
                        e              = kron(varargin{i}.Polytope,ones(sa1,1));
                        obj.Polytope   = [d,e];
                    end
            end
            if isempty(varargin{i}.UncertaintyMap)
                obj.UncertaintyMap     = {};
            else 
                if isempty(obj.UncertaintyMap)
                    obj.UncertaintyMap = varargin{i}.UncertaintyMap;
                else
                    [ma1,ma2]          = size(obj.UncertaintyMap{1});
                    [mb1,mb2]          = size(varargin{i}.UncertaintyMap{1});
                    la                 = length(obj.UncertaintyMap);
                    lb                 = length(varargin{i}.UncertaintyMap);
                    for ij = 1:la
                        H{ij}          = blkdiag(obj.UncertaintyMap{ij},zeros(mb1,mb2));
                    end
                    for ij = 1:lb
                        H{la+ij}       = blkdiag(zeros(ma1,ma2),varargin{i}.UncertaintyMap{ij});
                    end
                    obj.UncertaintyMap = H;
                end
            end
            a1 = length(varargin{i}.InputChannel);
            a2 = length(varargin{i}.OutputChannel);
            b1 = length(obj.InputChannel);
            b2 = length(obj.OutputChannel);
            if a1 ~= a2 || b1 ~= b2
                error('Error: Each delta-block should have the same number of in- and output channels (which may have a different dimension).');
            end
            if a1 > 1
                Q1  = aug(varargin{i}.LinNonlin,a1);
                Q2  = aug(varargin{i}.TimeInvTimeVar,a1);
                Q3  = aug(varargin{i}.StaticDynamic,a1);
                Q4  = aug(varargin{i}.Structure,a1);
                Q5  = aug(varargin{i}.Bounds,a1);
                Q6  = aug(varargin{i}.RateBounds,a1);
                Q7  = aug(varargin{i}.NormBounds,a1);
                Q8  = aug(varargin{i}.SectorBounds,a1);
                Q9  = aug(varargin{i}.SlopeBounds,a1);
                Q10 = aug(varargin{i}.DelayTime,a1);
                Q11 = aug(varargin{i}.DelayType,a1);
                Q12 = aug(varargin{i}.Passive,a1);
                Q13 = aug(varargin{i}.Odd,a1);
            else
                Q1  = varargin{i}.LinNonlin;
                Q2  = varargin{i}.TimeInvTimeVar;
                Q3  = varargin{i}.StaticDynamic;
                Q4  = varargin{i}.Structure;
                Q5  = varargin{i}.Bounds;
                Q6  = varargin{i}.RateBounds;
                Q7  = varargin{i}.NormBounds;
                Q8  = varargin{i}.SectorBounds;
                Q9  = varargin{i}.SlopeBounds;
                Q10 = varargin{i}.DelayTime;
                Q11 = varargin{i}.DelayType;
                Q12 = varargin{i}.Passive;
                Q13 = varargin{i}.Odd;
            end
            for j = 1:a1
                obj.InputChannel{b1+j}         = varargin{i}.InputChannel{j};
                obj.OutputChannel{b1+j}        = varargin{i}.OutputChannel{j};
                if strcmp(obj.ChannelClass{1},'U')
                    obj.LinNonlin{b1+j}        = Q1{j};
                    obj.TimeInvTimeVar{b1+j}   = Q2{j};
                    obj.StaticDynamic{b1+j}    = Q3{j};
                    obj.Structure{b1+j}        = Q4{j};
                    if isempty(varargin{i}.Bounds)
                        obj.Bounds{b1+j}       = [];
                    else
                        obj.Bounds{b1+j}       = Q5{j};
                    end
                    if isempty(varargin{i}.RateBounds)
                        obj.RateBounds{b1+j}   = [];
                    else
                        obj.RateBounds{b1+j}   = Q6{j};
                    end
                    if isempty(varargin{i}.NormBounds)
                        obj.NormBounds{b1+j}   = [];
                    else
                        obj.NormBounds{b1+j}   = Q7{j};
                    end
                    if isempty(varargin{i}.SectorBounds)
                        obj.SectorBounds{b1+j} = [];
                    else
                        obj.SectorBounds{b1+j} = Q8{j};
                    end
                    if isempty(varargin{i}.SlopeBounds)
                        obj.SlopeBounds{b1+j}  = [];
                    else
                        obj.SlopeBounds{b1+j}  = Q9{j};
                    end
                    if isempty(varargin{i}.DelayTime)
                        obj.DelayTime{b1+j}    = [];
                    else
                        obj.DelayTime{b1+j}    = Q10{j};
                    end
                    if isempty(varargin{i}.DelayType)
                        obj.DelayType{b1+j}    = [];
                    else
                        obj.DelayType{b1+j}    = Q11{j};
                    end
                    if isempty(varargin{i}.Passive)
                        obj.Passive{b1+j}      = [];
                    else
                        obj.Passive{b1+j}      = Q12{j};
                    end
                    if isempty(varargin{i}.Odd)
                        obj.Odd{b1+j}          = [];
                    else
                        obj.Odd{b1+j}          = Q13{j};
                    end
                end
            end
        end
        if strcmp(obj.ChannelClass{1},'U')
            obj.LinNonlin      = red(obj.LinNonlin,1);
            obj.TimeInvTimeVar = red(obj.TimeInvTimeVar,1);
            obj.StaticDynamic  = red(obj.StaticDynamic,1);
            obj.Structure      = red(obj.Structure,1);
            obj.Passive        = red(obj.Passive,1);
            obj.Odd            = red(obj.Odd,1);
            obj.Bounds         = red(obj.Bounds,2);
            obj.RateBounds     = red(obj.RateBounds,2);
            obj.NormBounds     = red(obj.NormBounds,2);
            obj.SectorBounds   = red(obj.SectorBounds,2);
            obj.SlopeBounds    = red(obj.SlopeBounds,2);
            obj.DelayTime      = red(obj.DelayTime,2);
            obj.DelayType      = red(obj.DelayType,2);
        end
    end
end
end
function B = red(A,type)
    switch type
        case 1
            Q{1} = A{1};
            l    = length(A);
            if l > 1
                for i = 1:l-1
                    if strcmp(A{i},A{i+1})
                        j(i) = 1;
                    else
                        j(i) = 0;
                    end
                end
                if sum(j) == l-1
                    B = Q;
                else
                    B = A;
                end
            end
        case 2
            Q{1} = A{1};
            l    = length(A);
            if l > 1
                for i = 1:l-1
                    if isempty(A{i}) && isempty(A{i+1})
                        j(i) = 1;
                    else
                        j(i) = 0;
                    end
                end
                if sum(j) == l-1
                    B = Q;
                else
                    B = A;
                end
            end   
    end
end
function A = aug(A,a1)
    l = length(A);
    p = A{1};
    if l == 1
        for ijk = 1:a1
            A{ijk} = p;
        end
    end
end
