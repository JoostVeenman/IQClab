classdef iqcinvar < iqcprob
    
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
% Date:        26-04-2022
% 
% -------------------------------------------------------------------------
%
% Description: Class definition to setup an IQC invariance analysis
%              problem. 
%
%              Note: This is a subclass of the class iqcprob
%
% Syntax:      prob = iqcinvar(varargin)
%
% Usage:       iqcinvar defines the properties for performing an IQC
%              invariance analysis through the command iqcinvariance.
%
%              The properties are provided as output of the IQC invariance
%              analysis through the function "iqcinvariance".
%
% Usage:       iqcprob defines an IQC-problem class, which supports the
%              parsers LMIlab as well as Yalmip.
%
%              For "prob = iqcinvar(varargin)", the varargin inputs come in
%              pairs and can be defined as: 
%
%              prob = iqcinvar('name','prop1','value1','prop2','value2',..)
%
%              Alternatively, the properties can be set by defining the
%              structure:
%
%              iqcinvarOpt.prop1 = value1
%              iqcinvarOpt.prop2 = value2
%                        ...              
%              iqcinvarOpt.propN = valueN
% 
%              prob = iqcinvar(iqcprobOpt)
%
%              Finally, the properties can be set and retrieved by using:
%
%                   set(prob,'propX','valueX') and get(prob,'propX')
%
%              The properties that can be specified are:
%            
%                1.) 'alp' specifies the L2-bound on the external
%                    disturbance w\in D_alp = {w\in L_2: ||w|| \leq0}
%                    (Default: 1)
%
%                2.) 'Tx0' specifies the non-zero elements of the initial
%                    state x0 (Default: [])
%
%              The properties that are obtained by an IQC-invariance
%              analysis are (depending on the scenario ('InvarCase')):
%
%                1.) 'alp' being the bound on the hyper ellipsoid Hinv:
%                    x(t)\in{x\in R^n: x^T Hinv x \leq alp^2}
%
%                2.) 'Hinv' defining the hyper ellipsoid region:
%                    x(t)\in{x\in R^n: x^T Hinv x \leq alp^2}
%
%                3.) 'PeakGain' being the peak gain on the performance
%                    channels
%
%
% -------------------------------------------------------------------------

properties
    % L2-bound on the disturbance input w
    alp       double {mustBeReal,mustBeFinite} = 1;
    
    % Tx0 selects the (non-zero) elements of initial state x0
    Tx0       double {mustBeReal,mustBeFinite} = [];
end

% internal properties
properties
    % Inverse of the hyper ellipsoidal variable H
    Hinv      double {mustBeReal,mustBeFinite} = [];
    
    % Peak gain bounds on performance channel w->z
    PeakGain  double {mustBeReal,mustBeFinite} = [];
    
    % Terminal cost constraint Z and scaling Tz
    Z                                          = [];
    Tz                                         = [];
end

methods
    function obj  = iqcinvar(varargin)
        if nargin == 1
            if ~isstruct(varargin{1})
                error('Error: The second input argument should be defined by a structure (see "help iqcinvar" for further details)');
            end
            if isfield(varargin{1},'alp')
                if isreal(varargin{1}.alp) && isscalar(varargin{1}.alp) && varargin{1}.alp > 0
                    obj.alp = varargin{1}.alp;
                else
                    error('Error: The option "alp" should be defined as a real positive scalar (Default = 1).');
                end
            end
            if isfield(varargin{1},'Tx0')
                if isreal(varargin{1}.Tx0) && ismatrix(varargin{1}.Tx0)
                    obj.alp = varargin{1}.Tx0;
                else
                    error('Error: The option "Tx0" should be defined as a vector that selects the non-zero elements of the initial condition x0 (Default = []).');
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
                    case 'alp'
                        if isreal(varargin{j(i)}) && isscalar(varargin{j(i)}) && varargin{j(i)} > 0
                            obj.alp = varargin{j(i)};
                        else
                            error('Error: The option "alp" should be defined as a real positive scalar (Default = 1).');
                        end
                    case 'Tx0'
                        if isreal(varargin{j(i)}) && ismatrix(varargin{j(i)})
                            obj.Tx0 = varargin{j(i)};
                        else
                            error('Error: The option "Tx0" should be defined as a vector that selects the non-zero elements of the initial condition x0 (Default = []).');
                        end
                end
            end
        end
        % Initialize the LMI problem
        clear yalmip
        obj.lmi = [];
        setlmis([]);
    end
    function iqcprob = fZ(iqcprob,Z,Tz)
    % ---------------------------------------------------------------------
    % Description: Augment the terminal cost matrix variable Z and scaling
    %              matrix Tz
    %
    % Syntax:      iqcprob = fZ(iqcprob,Z,Tz)
    %
    % Usage:       Specify the IQC problem "iqcprob" as input as well as
    %              the to-be-augmented matrix variable Z and Tz.
    % ---------------------------------------------------------------------
        Zprev         = get(iqcprob,'Z');
        if ~isobject(Zprev)
            if isempty(Zprev)
                Zprev = iqcvar(iqcprob,[0,0],'symmetric');
            end
        end
        set(iqcprob,'Z',blkdiag(Zprev,Z));
        
        Tzprev        = get(iqcprob,'Tz');
        set(iqcprob,'Tz',blkdiag(Tzprev,Tz)); 
    end
end
end