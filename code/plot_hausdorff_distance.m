function plot_hausdorff_distance(imf, curvature_of_head)
[fh, xh] = ksdensity(curvature_of_head);
P = [xh; fh]';
hds = zeros(1, size(imf, 2));

for i = 1:size(imf, 2)
    [fa, xia] = ksdensity(sum(imf(:, setdiff(1:size(imf, 2), i)), 2));
    Q = [xia; fa]';
    hds(i) = HausdorffDist(P, Q);
end

figure;
plot(1:size(imf, 2), hds, '-o');
xlabel('IMF component removed');
ylabel('Hausdorff Distance');
title('Hausdorff Distance for Different IMF Combinations');
end