main_init;
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
rhobulkK = 0.01;
rhobulkD = 0.5;
pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 1;
vD = 1;
N = 0.11;

%% Setup the first fixed point
stst.kind = 'stst';
stst.parameter = [rhobulkK, rhobulkD, pE, wE, pI, wI, vK, vD, N];
% This values are approximations obtained from dde23.
E1_st=4.84493267843259;
E2_st=1.15067126224841;
I1_st=0.329106801770075;
I2_st=0.263338159420767;
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

contpar=ind_rhobulkK;
namepar='rhobulkK';

%% Initialize the branch
stst_branch0 = SetupStst(funcs,'x',[E1_st; E2_st; I1_st; I2_st],'parameter',stst.parameter,...
    'contpar',contpar,'max_step',[contpar,0.005],'min_bound',...
    [contpar 0.02],'max_bound',[contpar 0.98],...
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

legend(ax1,'location','north');
xlabel(ax1,namepar);
ylabel(ax1,'E1');
title(ax1,sprintf(append('Bifurcation diagram for ',namepar)));

%% Continuation of Hopf bifurcation in two parameters
fprintf('----- Hopf branch -----\n');

parameter_bd={'max_bound',[ind_rhobulkK, 0.98; ind_rhobulkD, 0.98],...
    'min_bound',[ind_rhobulkK, 0.02; ind_rhobulkD, 0.02],...
    'max_step',[0,0.1; ind_rhobulkK,5e-3; ind_rhobulkD,5e-3]};

% Point numbers of the two hopf bifurcation points.
num_hopf = br_getflags(stst_branch_wbifs,'hopf');
% It seems that 'SetupHopf' sets hopf_branch0 to be a 2-par branch, so that
% 'br_contn' knows to continue it this way, rather than as a steady-state
% branch.
[hopf_branch0,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(1),...
    'contpar', [ind_rhobulkK,ind_rhobulkD],...
    'dir', ind_rhobulkK, 'step', -0.001, parameter_bd{:});

figure('Name', sprintf('N=%g', N), 'NumberTitle','off');clf;
ax2=gca;
title(ax2, sprintf('N=%g', N));
xlim(ax2,[0,1]);
ylim(ax2,[0,1]);

hopf_branch0=br_contn(funcs,hopf_branch0,700);

% figure('Name','Generating 2-par bifurcation diagram (2)','NumberTitle','off');clf;
% ax3=gca;
% title(ax3,'Generating 2-par bifurcation diagram (2)');
% [hopf_branch1,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(2),...
%     'contpar', [ind_rhobulkK,ind_rhobulkD],...
%     'dir', ind_rhobulkK, 'step', -0.001, parameter_bd{:});
% hopf_branch1=br_contn(funcs,hopf_branch1,500);
% xlim(ax3,[0,1]);
% ylim(ax3,[0,1]);