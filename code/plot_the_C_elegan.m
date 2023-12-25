function plot_the_C_elegan(mcd,start_frame,end_frame)

global pixel2um unit2um
pixel2um = 1.6835;
unit2um = 0.05;

root_folder_path = 'D:\Desktop';
folder_name = ['from_' num2str(start_frame) '_to_' num2str(end_frame)];
folder_path = fullfile(root_folder_path,folder_name);
create_folder(folder_path);

for i = start_frame:end_frame


    %% absolut reference frame
    centerline = convert_coordinates_and_add_stage_position(mcd(i).SegmentedCenterline, mcd(i).StagePosition);
    boundary_A = convert_coordinates_and_add_stage_position(mcd(i).BoundaryA, mcd(i).StagePosition);
    boundary_B = convert_coordinates_and_add_stage_position(mcd(i).BoundaryB, mcd(i).StagePosition);

    plot_3_curve(centerline,boundary_A,boundary_B);
    xlabel("x (mm)");
    ylabel("y (mm)");
    title("Absolute Reference Frame");

    file_name = sprintf("frame_%d_absolute.png",i);
    file_path = fullfile(folder_path,file_name);
    saveas(gcf,file_path)
    close;

    %% relative reference frame
    % centerline = get_pixel_coordinates_in_the_video(mcd(i).SegmentedCenterline);
    % boundary_A = get_pixel_coordinates_in_the_video(mcd(i).BoundaryA);
    % boundary_B = get_pixel_coordinates_in_the_video(mcd(i).BoundaryB);
    % 
    % plot_3_curve(centerline,boundary_A,boundary_B);
    % xlabel("x (pixel)");
    % ylabel("y (pixel)");
    % title("Relative Reference Frame");
    % 
    % file_name = sprintf("frame_%d_relative.png",i);
    % file_path = fullfile(folder_path,file_name);
    % saveas(gcf,file_path)
    % close;

end

end

function plot_3_curve(centerline,boundary_A,boundary_B)
figure;
axis equal
hold on
scatter(centerline(1,:),centerline(2,:),'black');
scatter(boundary_A(1,:),boundary_A(2,:),'red');
scatter(boundary_B(1,:),boundary_B(2,:),'blue');
legend("centerline","boundary A","boundary B")
end