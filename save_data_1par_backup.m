% This is original save_data_1par.m
function save_data_1par_backup(namepar, file_id, N, par_stst, y_stst, label_y, nunst_stst, par_list, freq_list, T_list, amp_list, max_list, min_list, par_bf_hopf, y_bf_hopf)

if namepar=='N'
    file_name_stationary = sprintf('./results/data/N/bf_stationary_%s.txt', label_y);
    file_name_periodic = sprintf('./results/data/N/bf_periodic_%s.txt', label_y);
    file_name_bfpt = sprintf('./results/data/N/bf_bfpt_%s.txt', label_y);
    f_bf = fopen(file_name_stationary, 'w');
    for i=1:length(par_stst)
        fprintf(f_bf, '%.9g %.9g %.9g\n', par_stst(i), y_stst(i), nunst_stst(i));
    end
    fclose(f_bf);
    f_bf = fopen(file_name_periodic, 'w');
    for i=1:length(par_list)
        fprintf(f_bf, '%.9g %.9g %.9g %.9g %.9g %.9g\n', par_list(i), freq_list(i), T_list(i), amp_list(i), max_list(i), min_list(i));
    end
    fclose(f_bf);
    f_bf = fopen(file_name_bfpt, 'w');
    for i=1:length(par_bf_hopf)
        fprintf(f_bf, '%.9g %.9g\n', par_bf_hopf(i), y_bf_hopf(i));
    end
    fclose(f_bf);
else
    % Each file_id corresponds to a different N.
    file_name_stationary = sprintf('./results/data/%s/bf_stationary_%d.txt', namepar, file_id);
    file_name_periodic = sprintf('./results/data/%s/bf_periodic_%d.txt', namepar, file_id);
    file_name_bfpt = sprintf('./results/data/%s/bf_bfpt_%d.txt', namepar, file_id);
    f_bf = fopen(file_name_stationary, 'w');
    fprintf(f_bf, '%g\n', N);
    for i=1:length(par_stst)
        fprintf(f_bf, '%.9g %.9g %.9g\n', par_stst(i), y_stst(i), nunst_stst(i));
    end
    fclose(f_bf);
    f_bf = fopen(file_name_periodic, 'w');
    fprintf(f_bf, '%g\n', N);
    for i=1:length(par_list)
        fprintf(f_bf, '%.9g %.9g %.9g %.9g %.9g %.9g\n', par_list(i), freq_list(i), T_list(i), amp_list(i), max_list(i), min_list(i));
    end
    fclose(f_bf); 
    f_bf = fopen(file_name_bfpt, 'w');
    for i=1:length(par_bf_hopf)
        fprintf(f_bf, '%.9g %.9g\n', par_bf_hopf(i), y_bf_hopf(i));
    end
    fclose(f_bf); 
end



end