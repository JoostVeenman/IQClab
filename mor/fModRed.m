function Go = fModRed(Gi,w)
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
% Date:        13-03-2020
% 
% -------------------------------------------------------------------------
%
% Description: Function is a use interface that allows to perform model
%              reduction techniques on a given state-space or transfer
%              function model.
%
% Syntax:      Go = fPostPros(Gi,w)
%
% Usage:       Go = fPostPros(Gi,w) will guide the user through various
%              options that allows to perform various model reduction
%              techniques on the given LTI state-space or transfer function
%              model Gi. Here w=[w1,w2,...,wN] is the frequency vector
%              (frequency band), which specifies the plotting range.
%
%              The following optional steps can be respectively executed:
%
%               1.) Perform an initial model order reduction by balanced
%                   truncation (see the function "fMOR" for further
%                   details).
%               2.) Remove "fast" poles/zeros (see the function "fRFZP" for
%                   further details).
%               3.) Smoothen a certain frequency band (see the function
%                   "fSmooth" for further details).
%               4.) Replace a complex zero or pole pair with a one real
%                   zero or pole (see the function "f2C21R" for further
%                   details).
%               5.) Perform a final model oder reduction by balanced
%                   truncation (see the function "fMOR" for further
%                   details). 
%
% -------------------------------------------------------------------------

% Preparations
Gi   = ss(Gi);
Go   = Gi;

h1 = get(0,'CurrentFigure');
if isempty(h1)
    h1 = 1;
else
    h1 = h1.Number;
end

%--------------------------------------------------------------------------
% 1.) Perform an initial model order reduction by balanced truncation
%--------------------------------------------------------------------------
r1  = input('Would you like perform a controller order reduction [y/n]: ','s');

if r1 ~= 'y' && r1 ~= 'n'
    r1 = 'n';
end

if r1 == 'y'
    r2 = 'y';
    while r2 == 'y'
        Gor = fMOR(Go,w,[],1e-9,'off');
        Gor = fMOR(Gor,w);
        h1  = gcf;
        h1  = h1.Number;
        r2 = input('Would you like to repeat the process? [y/n]: ','s');
        if r2 == 'y' || r2 == 'n'
            if r2 == 'y'
                h1  = gcf;
                h1  = h1.Number;
                close(h1)
            end
        else
            r2 = 'n';
        end
    end
    Go = Gor;
    close(h1)
end

%--------------------------------------------------------------------------
% 2.) Remove "fast" poles/zeros
%--------------------------------------------------------------------------
r1    = input('Would you like to remove fast poles/zeros [y/n]: ','s');
if r1 ~= 'y' && r1 ~= 'n'
    r1 = 'n';
end

if r1 == 'y'
    r2 = 'y';
    while r2 == 'y'
        fBodePZ(w,h1,Go,'b-','Original model');
        Gor = fRFZP(Go);
        close(h1)
        fBodePZ(w,h1,Gor,'g--','Reduced model',Go,'b-','Original model');
        r2 = input('Would you like to repeat the process? [y/n]: ','s');
        if r2 == 'y' || r2 == 'n'
            r2 = r2;
        else
            r2 = 'n';
        end
        close(h1)
    end
    Go = Gor;
end

%--------------------------------------------------------------------------
% 3.) Smoothen a certain frequency band
%--------------------------------------------------------------------------
r1    = input('Would you like to smoothen the controller [y/n]: ','s');
if r1 ~= 'y' && r1 ~= 'n'
    r1 = 'n';
end

r2 = 'y';
while r2 == 'y'
    if r1 == 'y'
        r2 = 'y';
        while r2 == 'y'
            fBodemagPZ(w,h1,Go,'b-','Original model');
            Gos = fSmooth(Go);
            close(h1)
            fBodemagPZ(w,h1,Gos,'g--','Reduced model',Go,'b-','Original model');
            r2 = input('Would you like to repeat the process? [y/n]: ','s');
            if r2 == 'y' || r2 == 'n'
                r2 = r2;
            else
                r2 = 'n';
            end
            close(h1)
        end
        Go = Gos;
    end
    r2 = input('Would you like to smoothen another entry of the controller? [y/n]: ','s');
    if r2 == 'y' || r2 == 'n'
        if r2 == 'y'
            r2 = 'y';
        end
    else
        r2 = 'n';
    end
end
%--------------------------------------------------------------------------
% 4.) Replace a complex zero or pole pair with a one real zero or pole
%--------------------------------------------------------------------------
r1    = input('Would you like to replace a complex zero or pole pair with a one real zero or pole respectively [y/n]: ','s');
if r1 ~= 'y' && r1 ~= 'n'
    r1 = 'n';
end

if r1 == 'y'
    r2 = 'y';
    while r2 == 'y'
        fBodemagPZ(w,h1,Go,'b-','Original model');
        Gos = f2C21R(Go);
        close(h1)
        fBodemagPZ(w,h1,Gos,'g--','Reduced model',Go,'b-','Original model');
        r2 = input('Would you like to repeat the process? [y/n]: ','s');
        if r2 == 'y' || r2 == 'n'
            r2 = r2;
        else
            r2 = 'n';
        end
        close(h1)
    end
    Go = Gos;
end

%--------------------------------------------------------------------------
% 5.) Perform a final model oder reduction by balanced truncation
%--------------------------------------------------------------------------
r1  = input('Would you like perform a final controller order reduction [y/n]: ','s');

if r1 ~= 'y' && r1 ~= 'n'
    r1 = 'n';
end

if r1 == 'y'
    r2 = 'y';
    while r2 == 'y'
        Gor = fMOR(Go,w,[],1e-9,'off');
        Gor = fMOR(Gor,w);
        r2 = input('Would you like to repeat the process? [y/n]: ','s');
        if r2 == 'y' || r2 == 'n'
            if r2 == 'y'
                h1  = gcf;
                h1  = h1.Number;
                close(h1)
            end
        else
            r2 = 'n';
        end
    end
    Go = Gor;
end
% close(h1)
end