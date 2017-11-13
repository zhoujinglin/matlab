close all; clear; clc;
%读取图像
image = imread('lena512.bmp');
I = image(:);
image = double(image);
[M, N] = size(image);
%统计图像的概率直方图
probability = zeros(1,256);
for i = 1:M
  for j = 1:N
    probability(image(i, j) + 1) = probability(image(i, j) + 1) + 1;
  end
end
probability = probability./(M*N);
x = 0:255;
bar(x,probability(x + 1),1); xlim([0, 255]);
%计算图像的熵(entropy)
entropy = 0;
for i = 1:256
  if probability(i) ~= 0
    entropy = entropy - probability(i)*log2(probability(i));
  end
end
%进行霍夫曼编码
dict = huffmandict(x,probability); %生成字典
enco = huffmanenco(I,dict); %编码
[p, q] = size(enco);
average_code_lenth = p/(M*N);
%计算压缩率(compression ratio)和冗余度(redundancy)
compression_ratio = 8/average_code_lenth;
redundancy = 1-(entropy/average_code_lenth);
disp(['熵 = ', num2str(entropy), '，压缩率 = ', num2str(compression_ratio), '，冗余度 = ', num2str(redundancy)])
