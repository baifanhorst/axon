main_init_more; % Here main_init_more rather than main_init should be used.
addpath('C:\Users\bai_f\Documents\research\Axon_growth\DelayedFeedback\mine_dde23');
%% Initial values of bifurcation parameters
% tau = N/nu/(1-rho_bulk)
% N==1, nu==1, rho_bulk==1/2, then we have the following initial values. 
alphaK = 0.8;
alphaD = 0.01;
pE = 6;
wE = 5;
pI = 6;
wI = 5;
vK = 1;
vD = 1;
N = 1;

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

contpar=ind_alphaD;
namepar='alphaD';

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

%% Extract frequency and amplitude of oscillation
fprintf('----- Extract frequencies -----\n');

JK = vK*J(rhobulkK);
%JD = vD*J(rhobulkD);
tauK = delay(rhobulkK, N, vK);
%tauD = delay(rhobulkD, N, vD);
%delays = [tauK, tauD];

tspan = [0, 200];
t_start = 100;
t_end = 200;



length_par = length(par_stst);
freq_list = [];
max_list = [];
min_list = [];
amp_list = [];
par_list = [];

for i=1:length_par
    if nunst_stst(i)>0
        alphaD = par_stst(i);
        gammaD = alpha_to_gamma3(alphaD, dD, wD, vD);
        rhobulkD = density_bulk(alphaD, gammaD);
        JD = vD*J(rhobulkD);
        tauD = delay(rhobulkD, N, vD);
        delays = [tauK, tauD];
        
        sol = dde23(@rhs_dde23, delays, @history, tspan);
        [freq, amp, T, amp_max, amp_min] = extract_freq_amp(sol, 1, t_start, t_end);
        if amp<1e-2
            freq = 0;
        end
        par_list(end+1) = alphaD;
        freq_list(end+1) = freq;
        max_list(end+1) = amp_max;
        min_list(end+1) = amp_min;
        amp_list(end+1) = amp;
    end
end

figure('Name','freq','NumberTitle','off')
plot(par_list, freq_list, 'b.');
xlabel(namepar);
ylabel('Frequency');

figure('Name','amp','NumberTitle','off')
plot(par_list, amp_list, 'b.');
xlabel(namepar);
ylabel('Amplitude');


figure(fig_bf);
ax1 = gca;
hold on;
plot(ax1,par_list,max_list,'k.', 'DisplayName','max');
plot(ax1,par_list,min_list,'r.', 'DisplayName','min');
legend(ax1,'location','north');