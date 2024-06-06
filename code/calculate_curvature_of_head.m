function curvature_of_head = calculate_curvature_of_head(curvature_of_centerline_all)

magnification = 100;

% Calculate the mean curvature of each set across a sliding window of 10 columns
for i = 1:5
    CH(:,i) = mean(curvature_of_centerline_all(:,i:i+9),2); % 1-15, head
end

curvature_of_head = mean(CH,2).*magnification;

end