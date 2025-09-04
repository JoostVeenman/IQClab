% -------------------------------------------------------------------------
%
% Copyright:   This is copyrighted material owned by Novantec B.V.
% Terms:       IQClab is available under a Creative Commons
%              (Attribution-NoDerivatives 4.0 International (CC BY-ND 4.0))
%              license: https://creativecommons.org/licenses/by-nd/4.0/
%              For further information please visit iqclab.eu
%
% Author:      J.Veenman
% Date:        16-04-2022
% 
% -------------------------------------------------------------------------
% Demo_013:    Regional analysis
%
%              This demo demonstrates how to perform a regional
%              (invariance) analysis for the performance case 'e2z'. For
%              external disturbance signal with ||w||<1, the example
%              computes the smallest ellipsoid for the performance output z.
% -------------------------------------------------------------------------
clc;close all;

% Specify inputs
clc
disp('------------------------------------------------------------------');
disp('Demo_013: A regional analysis example:');
disp('------------------------------------------------------------------');
disp(' ');
disp('This demo demonstrates how to perform a regional (invariance)     ');
disp('analysis for the performance case e2z. For external disturbance   ');
disp('signal with ||w||<1, the example computes the smallest ellipsoid  ');
disp('for the performance output z.');
disp(' ');
disp('------------------------------------------------------------------');

% Define the plant
A               = [-0.97,2.2,2.36,3.45;-0.21,-0.8,5.2,-0.35;-2.56,-4.97,-0.75,-9.75;-3.64,0.2,9.68,-0.64];
B1              = [-0.62;-0.7;-1.42;0];
B2              = [-0.1;-0.32;-0.84;0];
C1              = [0,-0.36,0.36,-0.57];
C2              = [1.5,-0.11,0,0.93;0.1,0,0,0];
D11             = -1.14;
D12             = -1.76;
D21             = [0;0];
D22             = [0;0];
M               = ss(A,[B1,B2],[C1;C2],[D11,D12;D21,D22]);
da              = 0.5*(5 + - 0.6);
dd              = 0.5*(5 - - 0.6);
d               = ureal('d',0);
de              = da+dd*d;
N               = lft(de,M);
[G,D]           = lftdata(N);

% Perform the IQC analysis
k               = [1,2,3,4];

clear H
for i  = 1:length(k)
    u           = iqcdelta('u','InputChannel',1,'OutputChannel',1,'Bounds',[-1,1]);
    u           = iqcassign(u,'ultis','Length',k(i),'RelaxationType','PR');
    p           = iqcdelta('p','InputChannel',2,'OutputChannel',2:3,'ChannelClass','P','PerfMetric','e2z');
    iv          = iqcinvariance(G,{u,p},'SolChk','on','eps',1e-6,'alp',1);    
    H{i}        = iv.H;
end

% Plot the ellipsoidal regions
clu             = [1;lft(-1,G)];
C               = clu.c(2:end,:);
W               = lyap(clu.a,clu.b*clu.b');
Wi              = W^-1;
Hp              = C*W*C';
[V,D]           = eig(Hp);
la1             = sqrt(D(1,1));
la2             = sqrt(D(2,2));
t               = linspace(0,2*pi,100);
rell            = V*[la1*cos(t);la2*sin(t)]; 
K               = null(Wi-C'/Hp*C);
F               = -clu.b'*Wi;
[A,B,C,D]       = ssdata(clu);
clc             = ss(-A+B*F,-B,C+D*F,D);
ti              = linspace(0,100,10000);
sc              = 2*randn(2,25);
IC              = zeros(size(K,1),length(sc));
for i = 1:length(sc)
   IC(:,i)      = sc(1,i)*K(:,1)+sc(2,i)*K(:,2);
end

figure(1);
plot(rell(1,:),rell(2,:),'linewidth',2,'color','k','LineWidth',1);hold on

sty             = {'r-','r--','r:','r-.'};
for q = 1:length(H)
    Hi          = H{q}^-1;
    [V,D]       = eig(Hi);
    la1         = sqrt(D(1,1));
    la2         = sqrt(D(2,2));
    rell        = V*[la1*cos(t);la2*sin(t)]; 
    figure(1)
    plot(rell(1,:),rell(2,:),sty{q},'LineWidth',2);hold on
end

% plot worst case trajectories
for j = 1:size(IC,2)
    xi          = IC(:,j);
    x0          = xi/sqrt(xi'*Wi*xi);
    [out,~]     = lsim(clc,zeros(1,length(ti)),ti,x0);
    si2         = out(:,1)';
    sig         = si2(end:-1:1);
    tig         = [ti linspace(ti(end)+(ti(end)-ti(end-1)),2*ti(end),10000)];
    [out,~]     = lsim(clu,[sig zeros(1,10000)],tig);
    y2          = -out(:,2:3);
    plot(y2(:,1),y2(:,2),'linewidth',1);
end
legend('Exact (tight) solution','IQCs with basis length 1','IQCs with basis length 2',...
    'IQCs with basis length 3','IQCs with basis length 4','Worst-case dist. with ||w||<1','Location','SE')
grid on
xlabel('e_1');
ylabel('e_2');
axis([-11,11,-0.9,0.9]);
set(gcf,'Position',[300,400,700,350]);
fCutFig(1,1);