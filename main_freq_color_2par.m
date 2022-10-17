%% Generating color representing frequencies
fprintf('----- Extract frequencies -----\n');

rhobulkK_list = linspace(0.02,0.98,97);
rhobulkD_list = linspace(0.02,0.98,97);
tspan = [0, 200];
t_start = 100;
t_end = 200;
freq_mat = zeros(length(rhobulkK_list), length(rhobulkD_list));
for i=1:length(rhobulkK_list)
    for j=1:length(rhobulkD_list)
        rhobulkK = rhobulkK_list(i)
        rhobulkD = rhobulkD_list(j)
        JK = vK*J(rhobulkK);
        JD = vD*J(rhobulkD);
        
        tauK = delay(rhobulkK, N, vK);
        tauD = delay(rhobulkD, N, vD);
        delays = [tauK, tauD];
        sol = dde23(@rhs_dde23, delays, @history, tspan);
        [freq, amp, T, amp_max, amp_min] = extract_freq_amp(sol, 1, t_start, t_end);
        if amp<1e-2
            freq = 0;
        end
        freq_mat(i,j) = freq;        
    end
end

[X,Y] = meshgrid(rhobulkK_list, rhobulkD_list);
figure(fig_2par);
ax = gca;
hold on;
fig = pcolor(X,Y,freq_mat);
set(fig, 'EdgeColor', 'none');
colorbar(ax);