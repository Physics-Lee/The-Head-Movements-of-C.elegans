function save_all_figures

% Ask user to select a directory for saving the figures
folder_name = uigetdir;

% Create the folder if it does not exist
create_folder(folder_name)

% Get handles to all open figures
fig_handles = findobj('Type', 'figure');

% Loop through all figure handles
for i = 1:length(fig_handles)
    figure(fig_handles(i)); % Make the figure active
    file_name_fig = sprintf('figure_%d.png', fig_handles(i).Number); % Create filename based on figure number
    full_file_path = fullfile(folder_name, file_name_fig); % Create full file path
    saveas(fig_handles(i), full_file_path); % Save figure to file
end

end