clear; close all; clc;

image = im2double(rgb2gray(imread('lena256.jpg')));
[M, N] = size(image);
sobel_x = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
sobel_y = [-1, -2, -1; 0, 0, 0; 1, 2, 1];

image_sobel_x = imfilter(image, sobel_x, 'conv', 0, 'same');
image_sobel_y = imfilter(image, sobel_y, 'conv', 0, 'same');
image_sobel = sqrt(image_sobel_x.^2 + image_sobel_y.^2);
figure(1);
subplot(2, 2, 1); imshow(image); title('原图像');
subplot(2, 2, 2); imshow(image_sobel_x); title('X方向梯度');
subplot(2, 2, 3); imshow(image_sobel_y); title('Y方向梯度');
subplot(2, 2, 4); imshow(image_sobel); title('边缘检测结果');

% level = graythresh(image);
% image_sobel_binary = imbinarize(image_sobel, level);
image_sobel_binary = zeros(M, N); level = 0.9;
for i = 1:M
  for j = 1:N
    if image_sobel(i, j) > level
      image_sobel_binary(i, j) = 1;
    else
      image_sobel_binary(i, j) = 0;
    end
  end
end

matlab_image_sobel = edge(image, 'sobel');
figure(2);
subplot(1, 2, 1); imshow(image_sobel_binary); title('自写程序边缘检测');
subplot(1, 2, 2); imshow(matlab_image_sobel); title('自带函数边缘检测');
