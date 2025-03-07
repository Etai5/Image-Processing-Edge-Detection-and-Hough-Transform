%% Image Processing Homework 4:
% Etai Wigman - 315875385      Itzik Nsimov - 206697849
disp("Etai Wigman 315875385     Itzik Nsimov 206697849")
%% 1) Edge Detection
clear variables
close all
clc
%% 1.1) Reading the Image
close all
clc
cm = imread("cameraman.tif");
cm = double(cm);
min = min(cm,[],"all");
max = max(cm,[],"all");
cm = (cm-min)/(max-min);
imshow(cm); title("Original Cameraman Image","FontSize",16);

clear max min

%% 1.2) Prewitt Edge Detector
% 1.2.2) Creating and Displaying two images with different thresholds
close all
clc
cm_p1 = dip_prewitt_edge(cm,0.2);
cm_p2 = dip_prewitt_edge(cm,0.05);

subplot(1,2,1)
imshow(cm_p1); title("Prewitt Edge Detector with Thresh=0.2",FontSize=16)
subplot(1,2,2)
imshow(cm_p2); title("Prewitt Edge Detector with Thresh=0.05",FontSize=16)

clear cm_p1 cm_p2

%% 1.3) Canny Edge Detector
% 1.3.2) using Canny with two different sets of parameters
close all
clc
cm_c1 = edge(cm,'canny');
cm_c2 = edge(cm,'canny',[0.1 0.4],2);

subplot(1,2,1)
imshow(cm_c1); title("Canny Edge Detector with Default Parameters",FontSize=16)
subplot(1,2,2)
imshow(cm_c2);
title("Canny Edge Detector with High/Low Threshold = 0.4/0.1 and \sigma=2",FontSize=16)

clear cm_c1 cm_c2

%% 1.3.3) Example of extraction of the thresh from default Canny
close all
clc
[~, thresh] = edge(cm,'canny');
% we can now use thresh to fine-tune our usage of Canny edge detector
edge(cm,'canny',thresh.*2); title("Canny with manipulation to found thresh",FontSize=16)

clear thresh

%% 2) Hough Transform
clear variables
close all
clc
%% 2.1 Hough line transform
close all
clc
%% a) Read and convert to normalized grayscale
floor = imread("floor.jpg");
floor = double(rgb2gray(floor));
min = min(floor,[],"all");
max = max(floor,[],"all");
floor = (floor-min)/(max-min);
imshow(floor); title("Normalized Grayscale Floor Image","FontSize",16);

clear max min

%% b) Extract Edges
close all
clc
BW = edge(floor);   % default is Sobel edge detection method
imshow(BW); title("Default edge function: Sobel",FontSize=16);

%% c+d) Hough Matrix for lines
close all
clc
[M1, Radius1, Theta1] = dip_hough_lines(BW,1,1);
[M2, Radius2, Theta2] = dip_hough_lines(BW,5,4);
subplot(1,2,1)
imshow(M1,[]); title("BW Hough Matrix using (R_0,\theta_0)=(1,1)",FontSize=16)
subplot(1,2,2)
imshow(M2,[]); title("BW Hough Matrix using (R_0,\theta_0)=(5,4)",FontSize=16)

%% e) Significant lines
close all
clc
% For (R0, Θ0) = (1, 1)
numpeaks = 4;
peaks1 = houghpeaks(M1, numpeaks);
imshow(floor);
hold on;
for k = 1:length(peaks1)
    rho = Radius1(peaks1(k, 1));
    theta = Theta1(peaks1(k, 2));
    x0 = rho * cosd(theta);
    y0 = rho * sind(theta);
    x1 = x0 + 1000 * (-sind(theta));
    y1 = y0 + 1000 * (cosd(theta));
    x2 = x0 - 1000 * (-sind(theta));
    y2 = y0 - 1000 * (cosd(theta));
    plot([x1, x2], [y1, y2], linewidth=2, color='r');
    title("4 most significant lines for (R_0,\theta_0)=(1,1)",FontSize=16);
end
hold off;

%%
% For (R0, Θ0) = (5, 4)
numpeaks = 4;
peaks2 = houghpeaks(M2, numpeaks,"Threshold",50);
imshow(floor);
hold on;
for k = 1:length(peaks2)
    rho = Radius2(peaks2(k, 1));
    theta = Theta2(peaks2(k, 2));
    x0 = rho * cosd(theta);
    y0 = rho * sind(theta);
    x1 = x0 + 1000 * (-sind(theta));
    y1 = y0 + 1000 * (cosd(theta));
    x2 = x0 - 1000 * (-sind(theta));
    y2 = y0 - 1000 * (cosd(theta));
    plot([x1, x2], [y1, y2], linewidth=2, color='r');
    title("4 most significant lines for (R_0,\theta_0)=(5,4)",FontSize=16);
end
hold off;

%% 2.2 Hough circle transform
clear variables
close all
clc
%% a) read coffee.png and turn it to grayscale and normalize it
close all
clc
coffee = imread("coffee.jpg");
coffee = double(rgb2gray(coffee));
min = min(coffee,[],"all");
max = max(coffee,[],"all");
coffee = (coffee-min)/(max-min);
imshow(coffee); title("Normalized Grayscale coffee Image","FontSize",16);

clear max min

%% b) Extract Edges
close all
clc
BW = edge(coffee);   % default is Sobel edge detection method
imshow(BW); title("Default edge of Coffee image",FontSize=16);

%% c) Hough Matrix for circles
close all
clc
[M1, Radius1] = dip_hough_circles(BW,1,1);
[M2, Radius2] = dip_hough_circles(BW,4,10);

%% d) Measuring Run-Time
close all
clc
tic; [original, Radius_original] = dip_hough_circles(BW,1,1); toc;

tic; [ref, Radius_ref] = dip_hough_circles(BW,2,16); toc;

%% e) Showing the Hough Matrix for circles
close all
clc
subplot(1,3,1)
imshow(M1(:,:,1),[]);
title("BW Hough Matrix for circles using (R_0,\theta_0)=(1,1)")
subplot(1,3,2)
imshow(M2(:,:,1),[]);
title("BW Hough Matrix for circles using (R_0,\theta_0)=(4,10)")

subplot(1,3,3)
imshow(ref(:,:,1),[]);
title("BW Hough Matrix for circles using (R_0,\theta_0)=(2,16)")

%% f) 5 most significant circles
close all
clc
numpeaks = 5;
% For (R0, Θ0) = (1, 1)
subplot(1,3,1)
peaks1 = houghpeaks3d(M1,numpeaks);
imshow(coffee);
hold on
viscircles(peaks1(:,[1,2]),Radius1(peaks1(:,3)));
title("5 most significant circles for (R_0,\theta_0)=(1,1)");
hold off;

% For reference (R0, Θ0) = (2, 16)
subplot(1,3,2)
peaks_ref = houghpeaks3d(ref,numpeaks);
imshow(coffee);
hold on
viscircles(peaks_ref(:,[1,2]),Radius_ref(peaks_ref(:,3)));
title("5 most significant circles for (R_0,\theta_0)=(2,16)");
hold off;

% For (R0, Θ0) = (4, 10)
subplot(1,3,3)
peaks2 = houghpeaks3d(M2,numpeaks);
imshow(coffee);
hold on
viscircles(peaks2(:,[1,2]),Radius2(peaks2(:,3)));
title("5 most significant circles for (R_0,\theta_0)=(4,10)");
hold off;