function plot_hausdorff_distance(imf, curvature_of_head)

% Density estimation of the curvature of the head
[fh, xh] = ksdensity(curvature_of_head);
P = [xh; fh]'; % Prepare the reference set P for distance calculation

% Initialize the Hausdorff distances array
hds = zeros(1, size(imf, 2));

% Calculate the Hausdorff distance by excluding each IMF one at a time
for i = 1:size(imf, 2)
    % Density estimation of the sum of all IMFs except the current one
    [fa, xia] = ksdensity(sum(imf(:, setdiff(1:size(imf, 2), i)), 2));
    Q = [xia; fa]'; % Prepare the set Q for distance calculation

    % Compute the Hausdorff distance between sets P and Q
    hds(i) = HausdorffDist(P, Q);
end

% Plot the Hausdorff distances
figure;
plot(1:size(imf, 2), hds, '-o');
xlabel('IMF component removed');
ylabel('Hausdorff Distance');
title('Hausdorff Distance for Different IMF Combinations');

end