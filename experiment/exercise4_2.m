clear; close all; clc;

path = 'f16takeoff_396s/';
video = zeros(240, 360, 396);
[X, Y, Z] = size(video);
for count = 1:Z % 读取图像，存入三维数组
  image_count = num2str(count);
  if count <= 9
    image_name = strcat('00', image_count, '.jpg');
  elseif count <= 99
    image_name = strcat('0', image_count, '.jpg');
  else
    image_name = strcat(image_count, '.jpg');
  end
  image = im2double(rgb2gray(imread([path image_name])));
  video(:, :, count) = image;
end
template_location_x = 90; template_location_y = 130; template_size_x = 40; template_size_y = 90;
search_x = 10; search_y = 10;
for k = 1:(Z - 1)
  absolute = zeros(2*search_x + 1, 2*search_y + 1);
  for x = -search_x:search_x
    for y = -search_y:search_y
      for i = template_location_x:(template_location_x + template_size_x)
        for j = template_location_y:(template_location_y + template_size_y)
          absolute(x + search_x + 1, y + search_y + 1) = absolute(x + search_x +1, y + search_y + 1) + abs(video(i, j, k) - video(i + x, j + y, k + 1));
        end
      end
    end
  end
  absolute_min = min(min(absolute));
  [x, y] = find(absolute == absolute_min);
  template_location_x = template_location_x + x - search_x - 1; template_location_y = template_location_y + y - search_y -1;
  % 输出视频
  imshow(video(:, :, k)); hold on; title(['Frame NO.',num2str(k)]);
  rectangle('Position', [template_location_y, template_location_x, template_size_y, template_size_x]);
  pause(0.01); hold off;
end
