function plot_psd_of_imfs(imf, time, curvature_of_centerline_all)
fs = size(curvature_of_centerline_all, 1) / (time(end) - time(1)); % Hz
figure;

for i = 1:size(imf, 2)
    [pxx, f] = periodogram(imf(:,i), rectwin(length(imf(:,i))), length(imf(:,i)), fs);
    subplot(4, 2, i);
    plot(f, pxx);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('PSD (W/Hz)');
    title(['Power Spectrum Density of imf ' num2str(i)]);
    xlim([0, max(instfreq(imf(:,i), fs))]);
end
end