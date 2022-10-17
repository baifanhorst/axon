%% Import biftool
clear;
close all; 
addpath('C:\dde_biftool\ddebiftool',...
    'C:\dde_biftool\ddebiftool_extra_psol',...
    'C:\dde_biftool\ddebiftool_extra_nmfm',...
    'C:\dde_biftool\ddebiftool_utilities');

%% Setup rhs
% There are two delays: tauK and tauD
funcs=set_funcs(...
    'sys_rhs', @rhs_more,...
    'sys_tau', @sys_tau_more,...   
    'sys_ntau',@()2);

%% Global parameters
% Global parameters will not be used as bifurcation parameters.
global pE dE1 dE2 wE nE KE betaE
global pI dI1 dI2 wI nI KI
global JK JD tauK tauD
global N
global vK vD
global wa_K wd_K wa_D wd_D K_K K_D R_K R_D
global M_K M_D eps
global dK wK dD wD
% global Vtip
global alphaK gammaK alphaD gammaD

%% Values of global parameters
dE1 = 1;
dE2 = 1;
dI1 = 1;
dI2 = 1;

nE = 4;
nI = 4;
KE = 2;
KI = 2;
betaE = 1;



wa_K = 5; wd_K = 1;
wa_D = 5; wd_D = 1;

K_K = wa_K/wd_K;
K_D = wa_D/wd_D;

R_K = 1/K_K;
R_D = 1/K_D;


M_K = 1;
M_D = 1;
eps = 1;

dK = 1;
wK = 1;
dD = 1;
wD = 1;


%% Indices of bifurcation parameters
ind_alphaK = 1;
ind_alphaD = 2;
ind_pE = 3;
ind_wE = 4;
ind_pI = 5;
ind_wI = 6;
ind_vK = 7;
ind_vD = 8;
ind_N = 9;

ind_dK = 10;
ind_wK = 11;
ind_dD = 12;
ind_wD = 13;

%% Functions for generating bifurcation diagrams
getpar=@(x,i)arrayfun(@(p)p.parameter(i),x.point);
getx=@(x,i)arrayfun(@(p)p.x(i),x.point);
% ???
bgetpar=@(x,i,bif)arrayfun(@(p)p.parameter(i),x.point(br_getflags(x,bif)));
bgetx=@(x,i,bif)arrayfun(@(p)p.x(i),x.point(br_getflags(x,bif)));
