clear; close all; clc;
% 读取红外图像
% path = '4a/';
% video_infrared = zeros(240, 320, 1506);
% [X_infrared, Y_infrared, T_infrared] = size(video_infrared);
% for count = 1 : 2 : T_infrared*2
%   image_count = num2str(count);
%   if count < 10
%     image_name = strcat('img_0000', image_count, '.bmp');
%   elseif count < 100
%     image_name = strcat('img_000', image_count, '.bmp');
%   elseif count < 1000
%     image_name = strcat('img_00', image_count, '.bmp');
%   else
%     image_name = strcat('img_0', image_count, '.bmp');
%   end
%   image = im2double(rgb2gray(imread([path image_name])));
%   video_infrared(:, :, (count + 1)/2) = image;
% end
% for i = 1 : T_infrared
%   imshow(video_infrared(:, :, i)); title(['Frame NO.', num2str(i)]); pause(0.01);
% end
% 读取彩色图像
path = '4b/';
video_colour = zeros(240, 320, 1507);
[X_colour, Y_colour, T_colour] = size(video_colour);
for count = 0 : 2 : (T_colour-1)*2
  image_count = num2str(count);
  if count < 10
    image_name = strcat('img_0000', image_count, '.bmp');
  elseif count < 100
    image_name = strcat('img_000', image_count, '.bmp');
  elseif count < 1000
    image_name = strcat('img_00', image_count, '.bmp');
  else
    image_name = strcat('img_0', image_count, '.bmp');
  end
  image = im2double(rgb2gray(imread([path image_name])));
  video_colour(:, :, (count + 2)/2) = image;
end
% for i = 1 : T_colour
%   imshow(video_colour(:, :, i)); title(['Frame NO.', num2str(i)]); pause(0.01);
% end
Sz = [240, 320, 1507];
% configuration , the number of direction at each size of pyramid
% total #direction is 3*Cf^2

% Cf = [3 3];

Cf = [6 12;6 6;6 6];

% default window parameter.
r = pi*[0.3 0.5 0.85 1.15];
alpha = 0.15;

% tic
% % 3 sec for size 64
% % 161 sec for size 256
% [F2, ind, cf] = ucurvwin3d_s(Sz, Cf, r, alpha);
% toc
%
% tic
% % 1 sec for size 64
% % 98 sec for size 256
% ydec = ucurvdec3d_s(video_colour, Cf, F2, ind, cf );
% toc

tic
ydec = ucurvdec3d(video_colour, Cf, r);
toc

tic
imr = ucurvrec3d(ydec, Cf, r);
toc
