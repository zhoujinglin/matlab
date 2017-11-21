clear; close all; clc;
% 读取图像
image = im2double(imread('homework4.jpg'));
[M, N] = size(image);
% 对图像进行最小值滤波
image_min = zeros(M, N);
s1 = zeros(1,9);
for i = 2:M-1
  for j = 2:N-1
    count = 0;
    for k = i - 1:i + 1
      for p = j - 1:j + 1
        count = count + 1;
        s1(count) = image(k, p);
      end
    end
    image_min(i, j) = min(s1);
  end
end
% 对图像二值化　取阈值为0.２
level = 0.2;
image_min_binary = zeros(M, N);
for i = 1:M
  for j = 1:N
    if image_min(i, j) > level
      image_min_binary(i, j) = 1;
    else
      image_min_binary(i, j) = 0;
    end
  end
end
subplot(1, 3, 1); imshow(image); title('原图像');
subplot(1, 3, 2); imshow(image_min); title('最小值滤波后的图像');
subplot(1, 3, 3); imshow(image_min_binary); title('二值化后的图像');
