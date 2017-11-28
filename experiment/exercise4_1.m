clear; close all; clc;

path = 'image200s/';
filecount = 200;
for i = 1:filecount
  name = num2str(i); % 读取图像，逐帧操作
  if i <= 9
    filename = strcat('0000000', name, '.bmp');
  elseif i <= 99
    filename = strcat('000000', name, '.bmp');
  elseif i <=199
    filename = strcat('00000', name, '.bmp');
  end
  image = im2double(imread([path filename]));
  % 将图像二值化，阈值为0.75
  [M, N] = size(image); level = 0.75;
  image_binary = zeros(M, N);
  for k = 1:M
    for j = 1:N
      if image(k, j) > level
        image_binary(k, j) = 1;
      else
        image_binary(k, j) = 0;
      end
    end
  end
  % 计算质心的位置
  size_image_binary = M*(N - 20) - sum(sum(image_binary(1:M, 1:(N - 20))));
  circle_x = 0; circle_y = 0;
  for k = 1:M
    for j = 1:(N - 20)
      circle_x = circle_x + k*(1 - image_binary(k, j));
      circle_y = circle_y + j*(1 - image_binary(k, j));
    end
  end
  circle_x = circle_x/size_image_binary; circle_y = circle_y/size_image_binary;
  % 输出目标跟踪视频
  imshow(image); hold on; plot (circle_y, circle_x, 'ro'); title(['Frame NO.' num2str(i)]);
  pause(0.01); hold off;
  % imshow([path filename]);
end
