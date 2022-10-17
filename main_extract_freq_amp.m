%% Extract frequency and amplitude of oscillation
fprintf('----- Extract frequencies -----\n');

length_par = length(par_stst);

par_list = [];
T_list = [];
freq_list = [];
amp_list = [];
max_list = [];
min_list = [];

amp_list_plot = [];
max_list_plot = [];
min_list_plot = [];



for i=1:5:length_par
    fprintf(append(namepar,sprintf('=%g\n',par_stst(i))));    
    if nunst_stst(i)>0
        switch namepar
            case 'pE'
                pE = par_stst(i);
            case 'pI'
                pI = par_stst(i);
            case 'wE'
                wE = par_stst(i);
            case 'wI'
                wI = par_stst(i);
            case 'rhobulkK'
                rhobulkK = par_stst(i);
            case 'rhobulkD'
                rhobulkD = par_stst(i);
            case 'vK'
                vK = par_stst(i);
            case 'vD'
                vD = par_stst(i);     
            case 'N'
                N = par_stst(i);
        end
     
        JK = vK*J(rhobulkK);
        JD = vD*J(rhobulkD);
        tauK = delay(rhobulkK, N, vK);
        tauD = delay(rhobulkD, N, vD);
        delays = [tauK, tauD];
        
        
        sol = dde23(@rhs_dde23, delays, @history, tspan);
        [freq, T, amp, amp_max, amp_min] = extract_freq_amp(sol, t_start, t_end, dt);
        
        par_list(end+1) = par_stst(i);
        freq_list(end+1) = freq;
        T_list(end+1) = T;
        % Note that amp_max, amp_min, amp are all lists.
        max_list = [max_list; amp_max];
        min_list = [min_list; amp_min];
        amp_list = [amp_list; amp];
        
        amp_list_plot(end+1) = amp(no_y);
        max_list_plot(end+1) = amp_max(no_y);
        min_list_plot(end+1) = amp_min(no_y);
        
        

    end
end

fig_freq=figure('Name','freq','NumberTitle','off');
plot(par_list, freq_list, 'b.');
xlabel(namepar);
ylabel('Frequency');
freq_max = max(freq_list);
freq_min = min(freq_list);
ylim([freq_min*0.8, freq_max*1.2]);

fig_period=figure('Name','Period','NumberTitle','off');
plot(par_list, T_list, 'b.');
xlabel(namepar);
ylabel('Period');
T_max = max(T_list);
T_min = min(T_list);
ylim([T_min*0.8, T_max*1.2]);






fig_amp=figure('Name','amp','NumberTitle','off');
plot(par_list, amp_list_plot, 'b.');
xlabel(namepar);
ylabel('Amplitude');


figure(fig_bf);
ax1 = gca;
hold on;
plot(ax1,par_list,max_list_plot,'k.', 'DisplayName','max');
plot(ax1,par_list,min_list_plot,'r.', 'DisplayName','min');
legend(ax1,'location','best');

