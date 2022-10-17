% This function saves the 2par branch in a txt file
% pointlist0 is 2par_branch0.point
% pointlist1 is 2par_branch1.point
% This is because the branch is extended in two directions.
function save_2par_branch(pointlist0, pointlist1, ind_x, ind_y, parname_x, parname_y)
    % Generate the file name
    filename = sprintf("./results/data_2par/2par_%s_%s_branch0.txt", parname_x, parname_y)
    f_bf = fopen(filename, 'w');
    n_point = length(pointlist0);
    for i = 1:n_point
        par_x = pointlist0(i).parameter(ind_x);
        par_y = pointlist0(i).parameter(ind_y);
        fprintf(f_bf, '%.9g %.9g\n', par_x, par_y);
    end
    fclose(f_bf);
    
    filename = sprintf("./results/data_2par/2par_%s_%s_branch1.txt", parname_x, parname_y)
    f_bf = fopen(filename, 'w');
    n_point = length(pointlist1);
    for i = 1:n_point
        par_x = pointlist1(i).parameter(ind_x);
        par_y = pointlist1(i).parameter(ind_y);
        fprintf(f_bf, '%.9g %.9g\n', par_x, par_y);
    end
    
    fclose(f_bf);
end