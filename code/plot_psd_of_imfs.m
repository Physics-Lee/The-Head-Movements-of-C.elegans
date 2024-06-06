function plot_psd_of_imfs(imf, time, curvature_of_centerline_all)

% Calculate the sampling frequency from the time array and curvature data
fs = size(curvature_of_centerline_all, 1) / (time(end) - time(1)); % Hz

% Create a new figure window
figure;

% Process each column in the 'imf' matrix, which corresponds to different IMFs
for i = 1:size(imf, 2)
    % Calculate the periodogram (PSD estimate) for the current IMF
    [pxx, f] = periodogram(imf(:,i), rectwin(length(imf(:,i))), length(imf(:,i)), fs);

    % Create a subplot for each IMF
    subplot(3, 2, i);
    plot(f, pxx);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('PSD (W/Hz)');
    title(['Power Spectrum Density of imf ' num2str(i)]);

    % Set x-axis limits based on the instantaneous frequency range of the current IMF
    xlim([0, max(instfreq(imf(:,i), fs))]);
end

end