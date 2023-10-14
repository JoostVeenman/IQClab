function W=fW(tmax,glf,type,Ts)
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
% Date:        20-08-2012
% 
% -------------------------------------------------------------------------
%
% Description: This function generates weighting filters that "cover" the
%              delay operator e^(-i\omega*tmax)-1, where tmax>0 is the
%              maximum delay.
%
% Syntax:      W = fW(tmax,glf,type)
%
% Usage:       Provide the inputs:
%              - "tmax": Maximum time-delay
%              - "glf" : Gain at low frequencies
%              - "type": There are three options:
%                 1.) Minimal gain W at high frequencies
%                 2.) Maximum bandwidth W
%                 3.) Both 1.) and 2.) at the cost of oder increase by 1.
%
%              If nargin > 3
%              - "Ts": Sample time. If Ts > 0, tmax should be an integer
%                                   denoting the number of delay steps
%
% -------------------------------------------------------------------------
% Example: Run the following code.
%
% tmax           = 0.02;
% om             = logspace(-4,3+log10(0.5/tmax),1000);
% Gdelay         = tf(1,1,'InputDelay',tmax)-1;
% glf            = 1e-3;
% theta          = 3.1608*pi/4;
% zer            = [-1.37/tmax -glf/tmax];
% pol            = (2*pi/(4*tmax))*[(exp(1i*theta)) (exp(1i*theta))'];
% P              = bodeoptions;
% P.PhaseVisible = 'off';
% P.FreqUnits    = 'Hz'; 
% P.Grid         = 'on';
% P.YLim         = [-70,20];
% s              = tf('s');
% W1             = 2*(s+0.5*glf/(1.5*tmax))/(s+1/(1.5*tmax));
% W2             = 4*(s+0.25*glf/(0.25*tmax))/(s+1/(0.25*tmax));
% W3             = zpk(zer,pol,2);
% bode(Gdelay,'g:',W1,'k-',W2,'b--',W3,'r-.',P,om);
% legend('Delay','Type 1','Type 2','Type 3','Location','NW')
% set(gcf,'Position',[23,44,1321,640]);
% -------------------------------------------------------------------------

if nargin < 4
    Ts = 0;
else
    if Ts > 0
        tmax = Ts*tmax;
    elseif Ts == -1
        error('The sampling time Ts must be a positive scalar');
    end
end

s=tf('s');
switch type
    case 1
        W     = 2*(s+glf/(3*tmax))/(s+2/(3*tmax));
    case 2
        W     = 4*(s+0.25*glf/(0.25*tmax))/(s+1/(0.25*tmax));
    case 3
        theta = 3.1608*pi/4;
        zer   = [-1.37/tmax -glf/tmax];
        pol   = (2*pi/(4*tmax))*[(exp(1i*theta)) (exp(1i*theta))'];
        W     = zpk(zer,pol,2);
end
if Ts > 0
    W         = c2d(W,Ts,'matched');
end
W             = balreal(ss(W));
end

