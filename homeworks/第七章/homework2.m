clc; clear; close all;

I_A = [-3, 0, 3, 3, 0; -4, 0, 4, 4, 0; -4, 0, 4, 4, 0; -4, 0, 4, 4, 0; -3, 0, 3, 3, 0];
[M, N] = size(I_A);
level = 2;
I_A_binary = zeros(M, N);
for i = 1:M
  for j = 1:N
    if I_A(i, j) > level
      I_A_binary(i, j) = 1;
    else
      I_A_binary(i, j) = 0;
    end
  end
end
subplot(1, 2, 1); imshow(I_A, []); title('原图像');
subplot(1, 2, 2); imshow(I_A_binary, []); title('二值化图像');
