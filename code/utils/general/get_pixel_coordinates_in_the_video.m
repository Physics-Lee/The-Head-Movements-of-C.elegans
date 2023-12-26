% Get pixel coordinates of a curve.
%
% Input: 200*1 numerical array in the form of xyxyxy...
% Output:2*100 numerical array where the 1st row is x and the 2nd row is y.
%
% 2023-12-25, Yixuan Li
%

function pixel_coordinates_in_the_video = get_pixel_coordinates_in_the_video(centerline)
pixel_coordinates_in_the_video = reshape(centerline, 2, 100);
end