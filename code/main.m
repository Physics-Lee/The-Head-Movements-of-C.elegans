dbstop if error;
clc;clear;close all;

%%
path = uigetdir;
list = get_all_files_of_a_certain_name_pattern_in_a_rootpath(path,"*mcd_corrected.mat");
[indx,tf] = listdlg('ListString',list,'ListSize',[800,600],'Name','Choose files');
for i = indx

    %% load your mcd
    full_path = list{i};
    mcd = load_data_from_mat(full_path);

    %% plot the C.elegan
    % start_frame = 19700;
    % end_frame = 19710;
    % plot_the_C_elegan(mcd,start_frame,end_frame);

    %% get the centerline in the absolut reference frame
    start_frame = 328;
    end_frame = 1243;
    [centerline_all,boundary_A_all,boundary_B_all] = get_centerlines_in_absolute_frame(mcd,start_frame,end_frame);

    %% get the curvature
    curvature_of_head = calculate_curvature(mcd,start_frame,end_frame);
end