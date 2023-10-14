close all;clc;

% user = 'new';
user = 'jv';

% IQCtoolbox paths
switch user
    case 'new'
        IQCtoolbox_path = 'xxx'; % Specify the IQC-toolbox path
        Yalmip_path     = 'yyy'; % Specify the yalmip toolbox path
        SDPT3_path      = 'zzz'; % Specify the SDPT3 solver path
%         Mosek_path      = 'nnn'; % Specify the Mosek solver path
    case 'jv'
        IQCtoolbox_path = 'C:\Data\IQClab';
        Yalmip_path     = 'C:\Data\yalmip';
        SDPT3_path      = 'C:\Data\sdpt3';
%         Mosek_path      = 'C:\Program Files\Mosek\9.1\toolbox\r2015a';
end

addpath(genpath(IQCtoolbox_path));
addpath(genpath(Yalmip_path));
addpath(genpath(SDPT3_path));
% addpath(genpath(Mosek_path));

clear user IQCtoolbox_path Yalmip_path SDPT3_path

run('IQClab_readme.m')