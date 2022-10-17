main_init;

%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 

pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 0.2;
vD = 1;

tauK = 1/vK/0.5;
tauD = 1/vD/0.5;
%% Setup the first fixed point
stst.kind = 'stst';
stst.parameter = [tauK, tauD, pE, wE, pI, wI, vK, vD];
% This values are approximations obtained from dde23.
E1_st=4.794;
E2_st=1.198;
I1_st=0.381;
I2_st=0.304;
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

contpar=ind_vK;

%% Initialize the branch
stst_branch0 = SetupStst(funcs,'x',[E1_st; E2_st; I1_st; I2_st],'parameter',stst.parameter,...
    'contpar',contpar,'max_step',[contpar,0.005],'min_bound',...
    [contpar 0.1],'max_bound',[contpar 15],...
    'newheuristics_tests',0);

%% Continuation of the branch
figure('Name','Generating branch for vK','NumberTitle','off');
clf;
ax1=gca;
title(ax1,sprintf('Generating branch for vK'));

[stst_branch0] = br_contn(funcs,stst_branch0,300);
%stst_branch0 = br_rvers(stst_branch0);
%[stst_branch0] = br_contn(funcs,stst_branch0,500);

[stst_branch_wbifs,stst_testfuncs]=LocateSpecialPoints(funcs,stst_branch0);
nunst_stst=GetStability(stst_branch_wbifs);

%% Plot bifurcation diagram
par_stst=getpar(stst_branch_wbifs,contpar);
E1_stst=getx(stst_branch_wbifs,1);
figure('Name','Bifurcation diagram for vK','NumberTitle','off');
clf;
ax1=gca;
plot(ax1,par_stst(nunst_stst==0),E1_stst(nunst_stst==0),'g.',...
    par_stst(nunst_stst==1),E1_stst(nunst_stst==1),'r.',...
    par_stst(nunst_stst==2),E1_stst(nunst_stst==2),'b.',...
    bgetpar(stst_branch_wbifs,contpar,'hopf'),bgetx(stst_branch_wbifs,1,'hopf'),'ks',...
    bgetpar(stst_branch_wbifs,contpar,'fold'),bgetx(stst_branch_wbifs,1,'fold'),'mo');
stst_lgtext={'unstable=0','unstable=1','unstable=2','hopf','fold'};
legend(ax1,stst_lgtext,'location','east');
xlabel(ax1,'vK');
ylabel(ax1,'E1');
title(ax1,sprintf('Bifurcation diagram for vK'));