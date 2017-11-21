clear; close all; clc;
%读取图像　取阈值为0.5
image = im2double(imread('starshape.jpg'));
[M, N] = size(image);
level = 0.5;
image_binary = imbinarize(image, level);
size_image_binary = M*N - sum(sum(image_binary));
% 计算质心的位置　并画出质心
circle_x = 0; circle_y = 0;
for i = 1:M
  for j = 1:N
    circle_x = circle_x + i*(1 - image_binary(i, j));
    circle_y = circle_y + j*(1 - image_binary(i, j));
  end
end
circle_x = circle_x/size_image_binary; circle_y = circle_y/size_image_binary;
subplot(1, 2, 1); imshow(image); hold on; plot(circle_y, circle_x, 'ro'); title('原图(含质心)');
subplot(1, 2, 2); imshow(image_binary); title('边缘检测');
