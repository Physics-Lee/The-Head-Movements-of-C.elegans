function plot_ch_minus_imf(curvature_of_head, imf)

figure;
subplot(3, 2, 1);
plot(curvature_of_head); % original signal
xlabel('frames');
ylabel('curvature*L');
title('curvature of head');

subplot(3, 2, 2);
plot(imf(:, 5) + imf(:, 4) + imf(:, 3) + imf(:, 2));
xlabel('frames');
ylabel('curvature*L');
title('IMF 5 + IMF 4 + IMF 3 + IMF 2');

subplot(3, 2, 3);
plot(imf(:, 5) + imf(:, 4) + imf(:, 3));
xlabel('frames');
ylabel('curvature*L');
title('IMF 5 + IMF 4 + IMF 3');

subplot(3, 2, 4);
plot(imf(:, 5) + imf(:, 4));
xlabel('frames');
ylabel('curvature*L');
title('IMF 5 + IMF 4');

subplot(3, 2, 5);
plot(imf(:, 5));
xlabel('frames');
ylabel('curvature*L');
title('IMF 5');

end