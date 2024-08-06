image_data = imread("IMG_20240406_110230.jpg");
image_data = rgb2gray(image_data);
timefilter = 1000;
bodyfilter = 1000;
magnification = 1;
h = fspecial('average', [timefilter bodyfilter]); % The algorithmic average value of the neighborhood around each pixel was calculated to smooth the data
image_data_filtered = imfilter(image_data * magnification,  h , 'replicate'); % N-D filtering of multidimensional images

figure;
imshow(image_data);

figure;
imshow(image_data_filtered);