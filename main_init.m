%% Import biftool

% [Clear] removes all variables from the current workspace, releasing them from system memory.
clear;
% [Close all] closes all figures whose handles are visible 
% Handle visibility is an advanced concept. Search HandleVisibility for more info.
% Currently, just think of [Close all] as closing all figures. 
close all; 

% Import BIFTOOL
% These paths are where BIFTOOL is saved.
addpath('C:\dde_biftool\ddebiftool',...
    'C:\dde_biftool\ddebiftool_extra_psol',...
    'C:\dde_biftool\ddebiftool_extra_nmfm',...
    'C:\dde_biftool\ddebiftool_utilities');

%% Setup rhs
% tauK and tauD will be the first two parameters, 
% as indicated by 'sys_tau'.

% For 'sys_ntau', we need to provide a function that specify the number of
% delays. In this case, there are two delays, so an anonymous function @() 2
% is used.
% For 'sys_tau', we need to provide a function that calculates the
% state-dependent delays. This function has a form of 
% sys_tau(delay_no, xx, par). [delay_no] specifies which delay is going to
% be calculated. In this case, delay_no = 1 or 2. [xx] is the state of the
% system, including E1(t), E2(t), I1(t), I2(t) as the 1st column,
% E1(t-tau1), E2(t-tau1), I1(t-tau1), I2(t-tau1) as the 2nd column, and 
% E1(t-tau2), E2(t-tau2), I1(t-tau2), I2(t-tau2) as the 3rd column.
% [par] is the parameter list. In this case,
% par = [rhobulkK, rhobulkD, pE, wE, pI, wI, vK, vD, N].
funcs=set_funcs(...
    'sys_rhs', @rhs,...
    'sys_tau', @sys_tau,...
    'sys_ntau',@()2);

%% Global parameters
% Global parameters will not be used as bifurcation parameters.
% pE, pI, wE, wI, JK, JD, tauK, tauD, N, vK, vD should not appear here.
% They appear here because the codes are copied from another script.
global pE dE1 dE2 wE nE KE betaE
global pI dI1 dI2 wI nI KI
global JK JD tauK tauD
global N
global vK vD
% These parameters are for the attachment and detachment. They will not be
% used.
% global wa_K wd_K wa_D wd_D K_K K_D R_K R_D
% global M_K M_D eps
global dK wK Vtip

%% Values of fixed parameters (compared to the bifurcation parameters)
dE1 = 1;
dE2 = 1;
dI1 = 1;
dI2 = 1;

nE = 4;
nI = 4;
KE = 2;
KI = 2;
betaE = 1;


% These parameters are for the attachment and detachment. They will not be
% used.
% wa_K = 5; wd_K = 1;
% wa_D = 5; wd_D = 1;
% 
% K_K = wa_K/wd_K;
% K_D = wa_D/wd_D;
% 
% R_K = 1/K_K;
% R_D = 1/K_D;


% M_K = 1;
% M_D = 1;
% eps = 1;

% dK = 1;
% Vtip = 1;

% wK = 1/Vtip;

%% Indices of bifurcation parameters
global ind_rhobulkK ind_rhobulkD ind_pE ind_wE ind_pI ind_wI
global ind_vK ind_vD ind_N

ind_rhobulkK = 1;
ind_rhobulkD = 2;
ind_pE = 3;
ind_wE = 4;
ind_pI = 5;
ind_wI = 6;
ind_vK = 7;
ind_vD = 8;
ind_N = 9;


%% Functions for generating bifurcation diagrams

% x represents a branch. x.point contains all the points along this branch.
% p represents a point.
% @(p)p.parameter(i) gets the ith parameter corresponding to point p.
% arrayfun applies the function @(p)... to all the points in x.point, and
% returns a list of values.
% Therefore, getpar gets the ith parameter of each point along the branch
% x.
% Similaryly, getx gets the ith coordinate of each point along the branch
% x.
% [getpar] and [getx] will be used when plotting bifurcation diagrams.
% [getpar] gets the parameter values and [getx] gets the y values. 
getpar=@(x,i)arrayfun(@(p)p.parameter(i),x.point);
getx=@(x,i)arrayfun(@(p)p.x(i),x.point);
% [br_getflags] returns the indices of bifucation points of type [bif].
% For example, br_getflags(stst_branch_wbifs, 'hopf') returns the indices
% of hopf bifucation points in the branch [stst_branch_wbifs]. Suppose that
% the result is [61, 157]. Then the 61st and 157th points will be selected
% by [x.point(br_getflags(x, bif))]. Then the anonymous function 
% [@(p) p.parameter(i)] returns the ith parameter for the two points.
% For example, if i=3, the values of pE at these two points will be
% returned. Note that for [bgetx], [i] represents different state
% variables. If i=3, then I1 will be returned.
bgetpar=@(x,i,bif)arrayfun(@(p)p.parameter(i),x.point(br_getflags(x,bif)));
bgetx=@(x,i,bif)arrayfun(@(p)p.x(i),x.point(br_getflags(x,bif)));
