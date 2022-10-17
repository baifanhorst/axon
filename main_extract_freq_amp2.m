%% Extract frequency and amplitude of oscillation
fprintf('----- Extract frequencies -----\n');

length_par = length(par_stst);
freq_list = [];
max_list = [];
min_list = [];
amp_list = [];
par_list = [];
T_list =[];

for i=1:length_par
    fprintf(append(namepar,sprintf('=%g\n',par_stst(i))));
    if nunst_stst(i)>0
        %[alphaK, alphaD, pE, wE, pI, wI, vK, vD, N, dK, wK, dD, wD]
        switch namepar
            case 'alphaK'
                alphaK = par_stst(i);
            case 'alphaD'
                alphaD = par_stst(i);
            case 'pE'
                pE = par_stst(i);
            case 'pI'
                pI = par_stst(i);
            case 'wE'
                wE = par_stst(i);
            case 'wI'
                wI = par_stst(i);
            case 'vK'
                vK = par_stst(i);
            case 'vD'
                vD = par_stst(i);     
            case 'N'
                N = par_stst(i);
            case 'dK'
                dK = par_stst(i);
            case 'wK'
                wK = par_stst(i);
            case 'dD'
                dD = par_stst(i);
            case 'wD'
                wD = par_stst(i);
        end
        
        gammaK = alpha_to_gamma3(alphaK, dK, wK, vK);
        rhobulkK = density_bulk(alphaK, gammaK);
        JK = vK*J(rhobulkK);
        tauK = delay(rhobulkK, N, vK);
        
        gammaD = alpha_to_gamma3(alphaD, dD, wD, vD);
        rhobulkD = density_bulk(alphaD, gammaD);
        JD = vD*J(rhobulkD);
        tauD = delay(rhobulkD, N, vD);
        
        
        delays = [tauK, tauD];
        
     
        JK = vK*J(rhobulkK);
        JD = vD*J(rhobulkD);
        tauK = delay(rhobulkK, N, vK);
        tauD = delay(rhobulkD, N, vD);
        delays = [tauK, tauD];
        
        sol = dde23(@rhs_dde23, delays, @history, tspan);
        [freq, amp, T, amp_max, amp_min] = extract_freq_amp(sol, 1, t_start, t_end, dt);
        if amp<1e-2
            freq = 0;
        end
        par_list(end+1) = par_stst(i);
        freq_list(end+1) = freq;
        max_list(end+1) = amp_max;
        min_list(end+1) = amp_min;
        amp_list(end+1) = amp;
        T_list(end+1) = T;
    end
end

fig_freq=figure('Name','freq','NumberTitle','off')
plot(par_list, freq_list, 'b.');
xlabel(namepar);
ylabel('Frequency');

fig_amp=figure('Name','amp','NumberTitle','off')
plot(par_list, amp_list, 'b.');
xlabel(namepar);
ylabel('Amplitude');

fig_period=figure('Name','Period','NumberTitle','off')
plot(par_list, T_list, 'b.');
xlabel(namepar);
ylabel('Period');


figure(fig_bf);
ax1 = gca;
hold on;
plot(ax1,par_list,max_list,'k.', 'DisplayName','max');
plot(ax1,par_list,min_list,'r.', 'DisplayName','min');
legend(ax1,'location','north');