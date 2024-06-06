function curvature_of_head = calculate_curvature_test(mcd, start_frame, end_frame)

%% Calculate curvatures of centerlines

n_frames = end_frame - start_frame + 1; % number of frames
n_curvpts = 100; % number of the points of the centerline

% Initialize variables
curvature_of_centerline_all = zeros(n_frames, n_curvpts);
time = zeros(n_frames, 1);

% Loop through frames
for j = 1:n_frames
    % Current frame
    i = start_frame + j - 1;

    % Get time
    time(j) = mcd(i).TimeElapsed;

    % Get centerline
    centerline = reshape(mcd(i).SegmentedCenterline, 2, []);

    % Calculate curvature
    curvature_of_centerline = calculate_the_curvature_of_a_centerline(centerline);

    % Save curvature data
    curvature_of_centerline_all(j, :) = curvature_of_centerline';
end

%% Plot the curvature diagram to verify the wave transmission (from head to tail)
plot_wave_transmission(curvature_of_centerline_all, time);

%% Calculate the curvature of head and body

% Calculate curvature
curvature_of_head = calculate_curvature_of_head(curvature_of_centerline_all);
curvature_of_body = calculate_curvature_of_body(curvature_of_centerline_all);

% Plot curvature of head and body
figure;
plot(time, curvature_of_head, 'red', time, curvature_of_body, 'blue');
xlabel('time (s)');
ylabel('curvature*L');
legend('curvature of head', 'curvature of body');

%% Variational Mode Decomposition (VMD) of the head

% Perform VMD
[imf, residual, info] = vmd(curvature_of_head);

% Plot each Intrinsic Mode Function (IMF)
figure;
num_imfs = size(imf, 2);
for i = 1:num_imfs
    subplot(3, 2, i);
    plot(time, imf(:, i));
    title(['IMF ' num2str(i)]);
    xlabel('time (s)');
    ylabel('curvature*L');
end

% Plot original signal
subplot(3, 2, 6);
plot(curvature_of_head, 'r');
xlabel('frames');
ylabel('curvature*L');
title('curvature of head');

% Combine high frequency IMFs
imf_high_f = imf(:, 2) + imf(:, 3) + imf(:, 4); % 1 is noise, so just 2+3+4
imf_low_f = imf(:, 5);

%% Plot Power Spectrum Density (PSD) of each IMF

% get sample frequency (Hz, or counts/s)
fs = size(curvature_of_centerline_all,1) / (time(end)-time(1)); % Hz

% 
figure;
num_imfs = size(imf, 2);
for i = 1:num_imfs
    [pxx, f] = periodogram(imf(:, i), rectwin(length(imf(:, i))), length(imf(:, i)), fs);
    subplot(3, 2, i);
    plot(f, 10*log10(pxx));
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB/Hz)');
    title(['Power Spectrum Density of IMF ' num2str(i)]);
end

end