function iqcprob = fBoundvar(iqcprob,var,bnd,scale)
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
% Date:        17-03-2020
%
% -------------------------------------------------------------------------
%
% Description: This function bounds the norm of an iqc variable "var" by
%              the constant or iqc variable "bnd". 
%
% Syntax:      iqcprob = fBoundvar(iqcprob,var,bnd)
%
% Usage:       iqcprob = fBoundvar(iqcprob,var,bnd) creates the following
%              LMI:
%
%              [0    , var] < bnd*I
%              [var^T, 0  ]
%
%              where bnd is either a positive scalar or an iqc variable.
%
%              This ensures ||var|| < bnd
%
% -------------------------------------------------------------------------
if nargin == 3
    if isobject(var) && isobject(bnd)
        n1      = var.Dim(1);
        n2      = var.Dim(2);
        m1      = bnd.Dim(1);
        m2      = bnd.Dim(2);
        if m1 ~= 1 || m2 ~= 1
           error('Error: The norm bound should be a defined as a scalar iqc variable');
        end   
        V       = blkdiag(oblkdiag(var),-bnd*eye(n1+n2));
        T       = [eye(n1+n2);eye(n1+n2)];
        iqcprob = iqclmi(iqcprob,V,-1,0,T);
    elseif isobject(var) && ~isobject(bnd)
        n1      = var.Dim(1);
        n2      = var.Dim(2);
        [m1,m2] = size(bnd);
        if m1 ~= 1 || m2 ~= 1
           error('Error: The norm bound should be defined as a positive scalar');
        end   
        V       = oblkdiag(var);
        T       = bnd*eye(n1+n2);
        iqcprob = iqclmi(iqcprob,V,-1,T);
    end
elseif nargin == 4
    if isobject(var) && isobject(bnd)
        n1      = var.Dim(1);
        n2      = var.Dim(2);
        m1      = bnd.Dim(1);
        m2      = bnd.Dim(2);
        if m1 ~= 1 || m2 ~= 1
           error('Error: The norm bound should be a defined as a scalar iqc variable');
        end   
        V       = blkdiag(oblkdiag(var),-bnd*eye(n1+n2));
        T       = [eye(n1+n2);sqrt(scale)*eye(n1+n2)];
        iqcprob = iqclmi(iqcprob,V,-1,0,T);
    elseif isobject(var) && ~isobject(bnd)
        n1      = var.Dim(1);
        n2      = var.Dim(2);
        [m1,m2] = size(bnd);
        if m1 ~= 1 || m2 ~= 1
           error('Error: The norm bound should be defined as a positive scalar');
        end   
        V       = oblkdiag(var);
        T       = scale*bnd*eye(n1+n2);
        iqcprob = iqclmi(iqcprob,V,-1,T);
    end
    
end 
end