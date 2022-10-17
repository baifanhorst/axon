%% Generating bifurcation diagram for density_bulk

%% Import biftool
clear;
close all; 
addpath('C:\dde_biftool\ddebiftool',...
    'C:\dde_biftool\ddebiftool_extra_psol',...
    'C:\dde_biftool\ddebiftool_extra_nmfm',...
    'C:\dde_biftool\ddebiftool_utilities');

%% Setup rhs
funcs=set_funcs(...
    'sys_rhs', @rhs,...
    'sys_tau', @()[1,2]);

%% Parameters except for tau
global dE1 dE2 nE KE betaE
%global pE dE1 dE2 wE nE KE betaE
global dI1 dI2 nI KI
% In this code, tauK, tauD, pE, wE, pI, wI, vK, vD are bifurcation parameters.
% They will be passed to the rhs.m, so they are not global.
% JK and JD depend on bulk densities (which depend on tauK and tauD) and vK
% and vD.
% Therefore, JK, JD are not global, either.

global N
%global vK vD
global wa_K wd_K wa_D wd_D K_K K_D R_K R_D
global M_K M_D eps
global dK wK Vtip

%% Set initial parameter values
% nE KE nI KI are set according to Bhargav's article
% All other parameters are set to be 1

pI = 6;
dE1 = 1;
dE2 = 1;
dI1 = 1;
dI2 = 1;

wI = 5;
nE = 4;
nI = 4;
KE = 2;
KI = 2;
betaE = 1;

N = 1;
vK = 1;
vD = 1;


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
Vtip = 1;

wK = vK/Vtip;




% alphaK = 0.3;
% gammaK = alpha_to_gamma(alphaK, v_K, M_K, eps, R_K);
% gammaK = alpha_to_gamma2(alphaK, dK, wK);
% alphaD = 0.7;
% gammaD = 0.1;

% density_bulk_K = density_bulk(alphaK, gammaK);
% density_bulk_K = 0.95;
% density_bulk_D = density_bulk(alphaD, gammaD);

% JK = vK*J(density_bulk_K);
% JD = vD*J(density_bulk_D);

%% Indices of parameters
ind_tauK = 1;
ind_tauD = 2;
ind_pE = 3;
ind_wE = 4;
ind_pI = 5;
ind_wI = 6;
ind_vK = 7;
ind_vD = 8;
%% Functions for generating bifurcation diagrams
getpar=@(x,i)arrayfun(@(p)p.parameter(i),x.point);
getx=@(x,i)arrayfun(@(p)p.x(i),x.point);
% ???
bgetpar=@(x,i,bif)arrayfun(@(p)p.parameter(i),x.point(br_getflags(x,bif)));
bgetx=@(x,i,bif)arrayfun(@(p)p.x(i),x.point(br_getflags(x,bif)));

%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
tauK = 1.052631578947368;
tauD = 2;
pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 1;
vD = 1;

%% Setup the first fixed point
stst.kind = 'stst';
stst.parameter = [tauK, tauD, pE, wE];
% This values are approximations obtained from dde23.
E1_st=4.84;
E2_st=1.15;
I1_st=0.33;
I2_st=0.26;
stst.x = [E1_st; E2_st; I1_st; I2_st];
% Correct the approximations above.
method = df_mthod(funcs, 'stst');
% If there is an error message caused by running the following line,
% you should check set_funcs to see whether the parameter positions are
% correct or not.
[stst, success] = p_correc(funcs, stst, [], [], method.point);
E1_st = stst.x(1);
E2_st = stst.x(2);
I1_st = stst.x(3);
I2_st = stst.x(4);

contpar=ind_tauK;

%% Initialize the branch for tauK
stst_branch0 = SetupStst(funcs,'x',[E1_st; E2_st; I1_st; I2_st],'parameter',stst.parameter,...
    'contpar',contpar,'max_step',[contpar,0.005],'min_bound',...
    [contpar N/vK],'max_bound',[contpar 100],...
    'newheuristics_tests',0);

%% Continuation of the branch for tauK
figure(1);
clf;
ax1=gca;
title(ax1,sprintf('steady states for tauK=%g',stst.parameter(ind_tauK)));
[stst_branch0] = br_contn(funcs,stst_branch0,2000);
[stst_branch_wbifs,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch0);
nunst_stst=GetStability(stst_branch_wbifs);

%% Plot bifurcation diagram for tauK
tauK_stst=getpar(stst_branch_wbifs,ind_tauK);
E1_stst=getx(stst_branch_wbifs,1);
figure(2);
clf;
ax1=gca;
plot(ax1,tauK_stst(nunst_stst==0),E1_stst(nunst_stst==0),'g.',...
    tauK_stst(nunst_stst==1),E1_stst(nunst_stst==1),'r.',...
    tauK_stst(nunst_stst==2),E1_stst(nunst_stst==2),'b.',...
    bgetpar(stst_branch_wbifs,ind_tauK,'hopf'),bgetx(stst_branch_wbifs,1,'hopf'),'ks',...
    bgetpar(stst_branch_wbifs,ind_tauK,'fold'),bgetx(stst_branch_wbifs,1,'fold'),'mo');
stst_lgtext={'unstable=0','unstable=1','unstable=2','hopf','fold'};
legend(ax1,stst_lgtext,'location','east');
xlabel(ax1,'tauK');
ylabel(ax1,'E1');
title(ax1,sprintf('steady states for tauK=%g',stst.parameter(ind_tauK)));
%% Plot bifurcation diagram for rho_bulk
rho_bulk_stst = 1 - N/vK./tauK_stst;
figure(3);
clf;
ax1=gca;
plot(ax1,rho_bulk_stst(nunst_stst==0),E1_stst(nunst_stst==0),'g.',...
    rho_bulk_stst(nunst_stst==1),E1_stst(nunst_stst==1),'r.',...
    rho_bulk_stst(nunst_stst==2),E1_stst(nunst_stst==2),'b.',...
    bgetpar(stst_branch_wbifs,ind_tauK,'hopf'),bgetx(stst_branch_wbifs,1,'hopf'),'ks',...
    bgetpar(stst_branch_wbifs,ind_tauK,'fold'),bgetx(stst_branch_wbifs,1,'fold'),'mo');
stst_lgtext={'unstable=0','unstable=1','unstable=2','hopf','fold'};
legend(ax1,stst_lgtext,'location','east');
xlabel(ax1,'bulk density');
xlim([0,1]);
ylabel(ax1,'E1');
title(ax1,sprintf('steady states for rho_bulk=%g',1-N/vK/stst.parameter(ind_tauK)));

%% Continuation of Hopf bifurcation in two parameters
fprintf('----- Hopf branch -----\n');

parameter_bd={'max_bound',[ind_tauK, 50; ind_tauD, 50],...
    'min_bound',[ind_tauK, N/vK; ind_tauD, N/vD],...
    'max_step',[0,0.1; ind_tauK,5e-3; ind_tauD,5e-3]};

% Point numbers of the two hopf bifurcation points.
num_hopf = br_getflags(stst_branch_wbifs,'hopf');
[hopf_branch0,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(1),...
    'contpar', [ind_tauK,ind_tauD],...
    'dir', ind_tauK, 'step', 0.002, parameter_bd{:});

figure(4);clf;
ax2=gca;
title(ax2,'Hopf in tauK-tauD plane');
hopf_branch0=br_contn(funcs,hopf_branch0,2000);

figure(5);clf;
ax3=gca;
title(ax3,'Hopf in tauK-tauD plane');
[hopf_branch1,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(2),...
    'contpar', [ind_tauK,ind_tauD],...
    'dir', ind_tauK, 'step', 0.002, parameter_bd{:});
hopf_branch1=br_contn(funcs,hopf_branch1,2000);

