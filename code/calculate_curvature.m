function curvature_of_head = calculate_curvature(mcd,start_frame,end_frame)

%% calculate curvatures of centerlines

n_frames = end_frame - start_frame + 1; % number of frames
n_curvpts = 100; % number of the points of the centerline

% init
curvature_of_centerline_all = zeros(n_frames,n_curvpts);
time = zeros(n_frames,1);

% loop
for j = 1:n_frames
    
    % frame now
    i = start_frame + j - 1;

    % get time
    time(j) = mcd(i).TimeElapsed;

    % get centerline
    centerline = reshape(mcd(i).SegmentedCenterline,2,[]);
    
    % calculate curvature
    curvature_of_centerline = calculate_the_curvature_of_a_centerline(centerline);

    % save
    curvature_of_centerline_all(j,:) = curvature_of_centerline';
    
end

%% plot the curvature diagram to verify the wave transmission (from head to tail)
plot_wave_transmission(curvature_of_centerline_all, time);

%% calculate the curvature of head and body

% calculate
curvature_of_head = calculate_curvature_of_head(curvature_of_centerline_all);
curvature_of_body = calculate_curvature_of_body(curvature_of_centerline_all);

% plot
figure;
plot(time,curvature_of_head,'red',time,curvature_of_body,'blue');
xlabel('time (s)');
ylabel('curvature*L');
legend('curvature of head','curvature of body');

%% vmd of the head

% vmd
[imf, ~, ~] = vmd(curvature_of_head); % Variational mode decomposition

% plot each imf
figure;
for i = 1:length(imf(1,:))    
    subplot(3,2,i);
    plot(time,imf(:,i));
    title(['IMF_' num2str(i)])
    xlabel('time (s)');
    ylabel('curvature*L');
end

% plot original signal
subplot(3,2,6)
plot(curvature_of_head, 'r');
xlabel('frames');
ylabel('curvature*L');
title('curvature of head')

% imf_high_f = imf(:,2) + imf(:,3) + imf(:,4); % 1 is noise, so just 2+3+4
imf_low_f = imf(:,5);

%% imf 5 vs curvature of the boday
figure;
plot(time, imf_low_f, 'red', time, curvature_of_body, 'blue');
xlabel('time (s)')
ylabel('curvature*L')
legend('IMF5','curvature of body')
title('IMF5 vs curvature of body')

%% plot original signal minus imfs
plot_ch_minus_imf(curvature_of_head, imf);

%% plot Power Spectrum Density of each imf
plot_psd_of_imfs(imf, time, curvature_of_centerline_all)

%% Hausdorff distance calculation and plot
plot_hausdorff_distance(imf, curvature_of_head)

%% compare ch and cb to know who is the low frequency part
compare_head_body(curvature_of_head, curvature_of_body);

%% save
save_all_figures;

end