function curvature_of_body = calculate_curvature_of_body(curvature_of_centerline_all)

magnification = 100;
curvature_of_body = mean(curvature_of_centerline_all(:,40:60),2).*magnification; % 40-60, body

end