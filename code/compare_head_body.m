function [lf_dynamics, hf_dynamics] = compare_head_body(head_dynamics, body_dynamics)

N = length(head_dynamics);
lf_template = zeros(N,3);

[imfv,~,~] = vmd(head_dynamics);

lf_template(:,1) = imfv(:,5);
lf_template(:,2) = imfv(:,4) + imfv(:,5);
lf_template(:,3) = imfv(:,3) + imfv(:,4) + imfv(:,5);
% candidate low frequency time sequence

head_dynamics_filtered = imfv(:,2) + imfv(:,3) + imfv(:,4) + imfv(:,5);
% noise is filtered

lf_template = lf_template - repmat(mean(lf_template,1),N,1);
body_dynamics = body_dynamics - mean(body_dynamics);
% substract mean

D = zeros(3,1);
figure;
for j = 1:3
    C = xcorr(body_dynamics, lf_template(:,j));
    [~,lag] = max(C(N:end));

    hd = lf_template(1:end-lag,j);
    bd = body_dynamics(1+lag:end);
    % shift the time sequence according the the cross correlation

    hd = hd/std(hd);
    bd = bd/std(bd);
    % scale the amplitude

    D(j) = sqrt(sum((hd - bd).^2)/(N-lag));
    % compute the mean square pairwise distance between points in the time sequence

    subplot(3,1,j);
    plot(1:N-lag,bd,'r','linewidth',1);
    hold on;
    plot(1:N-lag,hd, 'g', 'linewidth',1);
    box off
    xlabel("frame");
    ylabel("Normalized curvature");
    legend('body','head');
    switch j
        case 1
            title("IMF 5");
        case 2
            title("IMF 5 + IMF 4");
        case 3
            title("IMF 5 + IMF 4 + IMF 3");
    end
end

[~,I] = min(D);
lf_dynamics = lf_template(:,I);
hf_dynamics = head_dynamics_filtered - lf_dynamics;
switch I
    case 1
        disp('low frequency is imf5')
    case 2
        disp('low frequency is imf4 + imf5')
    case 3
        disp('low frequency is imf3 + imf4 + imf5')
end

end
