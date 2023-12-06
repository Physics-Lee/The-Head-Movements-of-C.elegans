clear;clc;close all;

% load human label file
option = 1;
switch option
    case 1
        human_label_file_name = 'human_label\human_label_N2.mat';
        start = 2;
        root_dir = 'F:\1_learning\research\taxis of C.elegans\data analysis of Colbert\data\Colbert\N2_taxis\NC20230310\w1';
        group = 'N2';
    case 2
        human_label_file_name = 'human_label\human_label_RIA.mat';
        start = 2;
        root_dir = 'F:\1_learning\research\taxis of C.elegans\data analysis of Colbert\data\Colbert';
        group = 'RIA-twk18';
end
load(human_label_file_name)
len = size(human_label, 1);
for i = start:len
    close all;
    
    % load mcd
    data_filename = 'mcd_corrected.mat';

    dir = fullfile(root_dir,data_filename);
    load(dir);
    
    % main
    mcd = eval(filename); % cool
    run_number = human_label{i,2};
    start_frame = human_label{i,3};
    end_frame = human_label{i,4};
    title_of_fig_1 = 'haha';
    pj_code_Yixuan_Li_modified(mcd,filename,title_of_fig_1,run_number,start_frame,end_frame)
    
end