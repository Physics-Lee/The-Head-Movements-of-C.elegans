dbstop if error;
clc;clear;close all;

%% super-parameter
global pixel2um unit2um
pixel2um = 1.683;
unit2um = 0.05;

%% choose folder
path = uigetdir;
list = get_all_files_of_a_certain_name_pattern_in_a_rootpath(path,"*mcd_corrected.mat");
[indx,tf] = listdlg('ListString',list,'ListSize',[800,600],'Name','Choose files');

%% you can use the longest forward in machine label or human label to plot
start_frame = 6200;
end_frame = 6800;

%% loop
for i = indx

    %% load your mcd
    full_path = list{i};
    mcd = load_data_from_mat(full_path);

    %% get the centerline in the absolut reference frame
    [centerline_all,boundary_A_all,boundary_B_all] = get_centerlines_in_absolute_frame(mcd,start_frame,end_frame);

    %% get the curvature
    curvature_of_head = calculate_curvature(mcd, start_frame, end_frame);
    
end