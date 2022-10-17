%% Continuation of Hopf bifurcation in two parameters
fprintf('----- Hopf branch -----\n');




parameter_bd={'max_bound',[contpar_2D_1, max_bound_2D_1; contpar_2D_2, max_bound_2D_2],...
    'min_bound',[contpar_2D_1, min_bound_2D_1; contpar_2D_2, min_bound_2D_2],...
    'max_step',[contpar_2D_1,max_step_2D_1; contpar_2D_2,max_step_2D_2]};


% parameter_bd={'max_bound',[ind_rhobulkK, 0.98; ind_rhobulkD, 0.98],...
%     'min_bound',[ind_rhobulkK, 0.02; ind_rhobulkD, 0.02],...
%     'max_step',[0,0.1; ind_rhobulkK,5e-3; ind_rhobulkD,5e-3]};

% Point numbers of the two hopf bifurcation points.
num_hopf = br_getflags(stst_branch_wbifs,'hopf');
% It seems that 'SetupHopf' sets hopf_branch0 to be a 2-par branch, so that
% 'br_contn' knows to continue it this way, rather than as a steady-state
% branch.
[hopf_branch0,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(1),...
    'contpar', [contpar_2D_1,contpar_2D_2],...
    'dir', contpar_2D_1, 'step', step_direction, parameter_bd{:});

fig_2par = figure('Name', figure_name, 'NumberTitle','off');clf;
ax2=gca;
hold on;
title(ax2, title_name);
xlim(ax2,x_range);
ylim(ax2,y_range);
xlabel(namepar_2D_1);
ylabel(namepar_2D_2);

hopf_branch0=br_contn(funcs,hopf_branch0,num_attempts);

% Continuation in the opposite direction
[hopf_branch1,suc] = SetupHopf(funcs, stst_branch_wbifs, num_hopf(1),...
    'contpar', [contpar_2D_1,contpar_2D_2],...
    'dir', contpar_2D_1, 'step', -step_direction, parameter_bd{:});
hopf_branch1=br_contn(funcs,hopf_branch1,num_attempts);