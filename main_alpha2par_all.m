main_init_more; % Here main_init_more rather than main_init should be used.
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
alphaK = 0.01;
alphaD = 0.2;
pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 1;
vD = 1;
N = 10;

dK = 1;
wK = 1;
dD = 1;
wD = 1;

%% Setup the first fixed point
stst.kind = 'stst';
stst.parameter = [alphaK, alphaD, pE, wE, pI, wI, vK, vD, N, dK, wK, dD, wD];
% This values are approximations obtained from dde23.

gammaK = alpha_to_gamma3(alphaK, dK, wK, vK);
rhobulkK = density_bulk(alphaK, gammaK);

gammaD = alpha_to_gamma3(alphaD, dD, wD, vD);
rhobulkD = density_bulk(alphaD, gammaD);


JK = vK*J(rhobulkK);
JD = vD*J(rhobulkD);
tauK = delay(rhobulkK, N, vK);
tauD = delay(rhobulkD, N, vD);
delays = [tauK, tauD];
tspan = [0, 200];
sol = dde23(@rhs_dde23, delays, @history, tspan);

E1_st=sol.y(1,end);
E2_st=sol.y(2,end);
I1_st=sol.y(3,end);
I2_st=sol.y(4,end);
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

contpar=ind_alphaK;
namepar='alphaK';

%% Initialize the branch
stst_branch0 = SetupStst(funcs,'x',[E1_st; E2_st; I1_st; I2_st],'parameter',stst.parameter,...
    'contpar',contpar,'max_step',[contpar,0.01],'min_bound',...
    [contpar 0.01],'max_bound',[contpar 1],...
    'newheuristics_tests',0);

%% Continuation of the branch
figure('Name',append('Generating branch for ',namepar),'NumberTitle','off');
clf;
ax1=gca;
title(ax1,sprintf(append('Generating branch for ',namepar)));
[stst_branch0] = br_contn(funcs,stst_branch0,500);
[stst_branch_wbifs,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch0);
nunst_stst=GetStability(stst_branch_wbifs);

%% Plot bifurcation diagram
par_stst=getpar(stst_branch_wbifs,contpar);
E1_stst=getx(stst_branch_wbifs,1);
fig_bf = figure('Name',append('Bifurcation diagram for ',namepar),'NumberTitle','off');
clf;
ax1=gca;
hold on;
plot(ax1,par_stst(nunst_stst==0),E1_stst(nunst_stst==0),'g.','DisplayName','stable');
plot(ax1,par_stst(nunst_stst>0),E1_stst(nunst_stst>0),'b.','DisplayName','unstable');
plot(ax1,bgetpar(stst_branch_wbifs,contpar,'hopf'),bgetx(stst_branch_wbifs,1,'hopf'),'ks','DisplayName','hopf');

legend(ax1,'location','best');
xlabel(ax1,namepar);
ylabel(ax1,'E1');
title(ax1,sprintf(append('Bifurcation diagram for ',namepar)));

%% Continuation of Hopf bifurcation in two parameters
fprintf('----- Hopf branch -----\n');

ind_par1 = ind_alphaK;
ind_par2 = ind_alphaD;
namepar1 = 'alphaK';
namepar2 = 'alphaD';

parameter_bd={'max_bound',[ind_par1, 1; ind_par2, 1],...
    'min_bound',[ind_par1, 0.01; ind_par2, 0.01],...
    'max_step',[ind_par1,0.01; ind_par2,0.01]};

% Point numbers of the two hopf bifurcation points.
num_hopf = br_getflags(stst_branch_wbifs,'hopf');
% It seems that 'SetupHopf' sets hopf_branch0 to be a 2-par branch, so that
% 'br_contn' knows to continue it this way, rather than as a steady-state
% branch.
[hopf_branch0,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(1),...
    'contpar', [ind_par1,ind_par2],...
    'dir', ind_par1, 'step', 0.01, parameter_bd{:});

figure('Name', sprintf('N=%g', N), 'NumberTitle','off');clf;
ax2=gca;
title(ax2, sprintf('N=%g', N));
xlim(ax2,[0,1]);
ylim(ax2,[0,1]);

hopf_branch0=br_contn(funcs,hopf_branch0,500);

hopf_branch0=br_rvers(hopf_branch0);
hopf_branch0=br_contn(funcs,hopf_branch0,500);

xlabel(ax2, namepar1);
ylabel(ax2, namepar2);

