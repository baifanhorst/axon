main_init;
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
rhobulkK = 0.5;
rhobulkD = 0.5;
pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 1;
vD = 1;
N = 0.2;

%% Setup the first fixed point
stst.kind = 'stst';
stst.parameter = [rhobulkK, rhobulkD, pE, wE, pI, wI, vK, vD, N];
% This values are approximations obtained from dde23.
E1_st=1.65103923106191;
E2_st=2.06386924257624;
I1_st=1.77117764962314;
I2_st=1.41696438084568;
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
    [contpar 0.01],'max_bound',[contpar 0.99],...
    'newheuristics_tests',0);

%% Continuation of the branch
figure('Name',append('Generating branch for ',namepar),'NumberTitle','off');
clf;
ax1=gca;
title(ax1,sprintf(append('Generating branch for ',namepar)));
[stst_branch0] = br_contn(funcs,stst_branch0,500);
stst_branch0 = br_rvers(stst_branch0);
[stst_branch0] = br_contn(funcs,stst_branch0,500);
[stst_branch_wbifs,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch0);
nunst_stst=GetStability(stst_branch_wbifs);

%% Plot bifurcation diagram
par_stst=getpar(stst_branch_wbifs,contpar);
E1_stst=getx(stst_branch_wbifs,1);
figure('Name',append('Bifurcation diagram for ',namepar),'NumberTitle','off');
clf;
ax1=gca;
plot(ax1,par_stst(nunst_stst==0),E1_stst(nunst_stst==0),'g.',...
    par_stst(nunst_stst==1),E1_stst(nunst_stst==1),'r.',...
    par_stst(nunst_stst==2),E1_stst(nunst_stst==2),'b.',...
    bgetpar(stst_branch_wbifs,ind_rhobulkK,'hopf'),bgetx(stst_branch_wbifs,1,'hopf'),'ks',...
    bgetpar(stst_branch_wbifs,ind_rhobulkK,'fold'),bgetx(stst_branch_wbifs,1,'fold'),'mo');
stst_lgtext={'unstable=0','unstable=1','unstable=2','hopf','fold'};
legend(ax1,stst_lgtext,'location','north');
xlabel(ax1,namepar);
ylabel(ax1,'E1');
title(ax1,sprintf(append('Bifurcation diagram for ',namepar)));


%% Continuation of Hopf bifurcation in two parameters
fprintf('----- Hopf branch -----\n');

parameter_bd={'max_bound',[ind_rhobulkK, 0.99; ind_rhobulkD, 0.99],...
    'min_bound',[ind_rhobulkK, 0.01; ind_rhobulkD, 0.01],...
    'max_step',[0,0.1; ind_rhobulkK,5e-3; ind_rhobulkD,5e-3]};

% Point numbers of the two hopf bifurcation points.
num_hopf = br_getflags(stst_branch_wbifs,'hopf');
% It seems that 'SetupHopf' sets hopf_branch0 to be a 2-par branch, so that
% 'br_contn' knows to continue it this way, rather than as a steady-state
% branch.
[hopf_branch0,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(1),...
    'contpar', [ind_rhobulkK,ind_rhobulkD],...
    'dir', ind_rhobulkK, 'step', -0.001, parameter_bd{:});

figure('Name','Generating 2-par bifurcation diagram (1)','NumberTitle','off');clf;
ax2=gca;
title(ax2,'Generating 2-par bifurcation diagram (1)');
hopf_branch0=br_contn(funcs,hopf_branch0,500);

figure('Name','Generating 2-par bifurcation diagram (2)','NumberTitle','off');clf;
ax3=gca;
title(ax3,'Generating 2-par bifurcation diagram (2)');
[hopf_branch1,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(2),...
    'contpar', [ind_rhobulkK,ind_rhobulkD],...
    'dir', ind_rhobulkK, 'step', -0.001, parameter_bd{:});
hopf_branch1=br_contn(funcs,hopf_branch1,500);