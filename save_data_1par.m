function save_data_1par(namepar, contpar, file_id, N, branchpoints, nunst_stst, par_list, freq_list, T_list, amp_list, max_list, min_list)




% branchpoints = stst_branch_wbifs.point
if namepar=='N'
    file_name_stationary = sprintf('./results/data/N/bf_stationary.txt');
    file_name_periodic = sprintf('./results/data/N/bf_periodic.txt');
    file_name_bfpt = sprintf('./results/data/N/bf_bfpt.txt');
    % Save the stationary branch
    f_bf = fopen(file_name_stationary, 'w');
    for i=1:length(branchpoints)
        fprintf(f_bf, '%.9g %.9g %.9g %.9g %.9g %.9g\n', branchpoints(i).parameter(contpar), branchpoints(i).x(1), branchpoints(i).x(2), branchpoints(i).x(3), branchpoints(i).x(4), nunst_stst(i));
    end
    fclose(f_bf);
    % Save the periodic branch
    f_bf = fopen(file_name_periodic, 'w');
    for i=1:length(par_list)
        fprintf(f_bf, '%.9g %.9g %.9g ', par_list(i), freq_list(i), T_list(i));
        fprintf(f_bf, '%.9g %.9g %.9g %.9g ', amp_list(i,1), amp_list(i,2), amp_list(i,3), amp_list(i,4));
        fprintf(f_bf, '%.9g %.9g %.9g %.9g ', max_list(i,1), max_list(i,2), max_list(i,3), max_list(i,4));
        fprintf(f_bf, '%.9g %.9g %.9g %.9g', min_list(i,1), min_list(i,2), min_list(i,3), min_list(i,4));
        fprintf(f_bf, '\n');
    end
    fclose(f_bf);
    % Save the bifurcation points
    f_bf = fopen(file_name_bfpt, 'w');
    for i=1:length(branchpoints)
        if ~isempty(branchpoints(i).flag) && strcmp(branchpoints(i).flag,'hopf')
            fprintf(f_bf, '%.9g %.9g %.9g %.9g %.9g\n', branchpoints(i).parameter(contpar), branchpoints(i).x(1), branchpoints(i).x(2), branchpoints(i).x(3), branchpoints(i).x(4));
        end
    end
    fclose(f_bf);
    
else
    % Each file_id corresponds to a different N.
    file_name_stationary = sprintf('./results/data/%s/bf_stationary_%d.txt', namepar, file_id);
    file_name_periodic = sprintf('./results/data/%s/bf_periodic_%d.txt', namepar, file_id);
    file_name_bfpt = sprintf('./results/data/%s/bf_bfpt_%d.txt', namepar, file_id);
    
    f_bf = fopen(file_name_stationary, 'w');
    fprintf(f_bf, '%g\n', N);
    for i=1:length(branchpoints)
        fprintf(f_bf, '%.9g %.9g %.9g %.9g %.9g %.9g\n', branchpoints(i).parameter(contpar), branchpoints(i).x(1), branchpoints(i).x(2), branchpoints(i).x(3), branchpoints(i).x(4), nunst_stst(i));
    end
    fclose(f_bf);
    
    f_bf = fopen(file_name_periodic, 'w');
    fprintf(f_bf, '%g\n', N);
    for i=1:length(par_list)
        fprintf(f_bf, '%.9g %.9g %.9g ', par_list(i), freq_list(i), T_list(i));
        fprintf(f_bf, '%.9g %.9g %.9g %.9g ', amp_list(i,1), amp_list(i,2), amp_list(i,3), amp_list(i,4));
        fprintf(f_bf, '%.9g %.9g %.9g %.9g ', max_list(i,1), max_list(i,2), max_list(i,3), max_list(i,4));
        fprintf(f_bf, '%.9g %.9g %.9g %.9g', min_list(i,1), min_list(i,2), min_list(i,3), min_list(i,4));
        fprintf(f_bf, '\n');
    end
    fclose(f_bf); 
    
    f_bf = fopen(file_name_bfpt, 'w');
    fprintf(f_bf, '%g\n', N);
    for i=1:length(branchpoints)
        if ~isempty(branchpoints(i).flag) && strcmp(branchpoints(i).flag,'hopf')
            fprintf(f_bf, '%.9g %.9g %.9g %.9g %.9g\n', branchpoints(i).parameter(contpar), branchpoints(i).x(1), branchpoints(i).x(2), branchpoints(i).x(3), branchpoints(i).x(4));
        end
    end
    fclose(f_bf); 
end



end