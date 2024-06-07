function plot_wave_transmission(curvature_of_centerline_all, time)

% smooth
timefilter = 5; % filter along y axis
bodyfilter = 10; % filter along x axis
magnification = 1;
h = fspecial('average', [timefilter bodyfilter]); % The algorithmic average value of the neighborhood around each pixel was calculated to smooth the data
curvdatafiltered = imfilter(curvature_of_centerline_all*magnification,  h , 'replicate'); % N-D filtering of multidimensional images

% plot
figure;
imagesc(curvdatafiltered); % imagesc is MATLAB function

% color
cmap = redgreencmap;
cmap(:,3)=cmap(:,2);
cmap(:,2)=0;
colormap(cmap);
colorbar;
clim([-0.1 0.1]);

% decorate
hold on;
title('cuvature diagram');
set(gca,'XTICK',[1 20 40 60 80 100]);
set(gca,'XTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);
y_tick=get(gca,'YTICK');
set(gca,'YTICKLABEL',time(y_tick));
xlabel('fractional distance along the centerline (head=0; tail=1)');
ylabel('time (s)');

end