clc; clear; close all;

A = [0, 1, 2, 3, 0, 1; 1, 2, 3, 0, 1, 2; 2, 3, 0, 1, 2, 3; 3, 0, 1, 2, 3, 0; 0, 1, 2, 3, 0, 1; 1, 2, 3, 0, 1, 2];
[M, N] = size(A);
L = 4;
G1 = zeros(L, L); G2 = zeros(L, L);
%　灰度共生矩阵G1
for x = 1:M - 1
  for y = 2:N
    G1(A(x, y) + 1, A(x + 1, y - 1) + 1) = G1(A(x, y) + 1, A(x + 1, y - 1) + 1) + 1;
  end
end
%　灰度共生矩阵G2
for x = 1:M - 1
  for y = 1:N - 1
    G2(A(x, y) + 1, A(x + 1, y + 1) + 1) = G2(A(x, y) + 1, A(x + 1, y + 1) + 1) + 1;
  end
end
