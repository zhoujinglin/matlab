clear; close all; clc;

x = 1024; y = 1024;
image = zeros(x, y);

disp('Take curvelet transform: fdct_wrapping');
tic; C = fdct_wrapping(image, 0, 2, 8, 64);

s = 7;
w = 10;
[A, B] = size(C{s}{w});
a = ceil((A + 1)/2);
b = ceil((B + 1)/5);
C{s}{w}(a, b) = 1;

disp('Take adjoint curvelet transform: ifdct_wrapping');
tic; Y = ifdct_wrapping(C, 0);

F = ifftshift(fft2(fftshift(Y)));
subplot(1, 2, 1); colormap gray; imagesc(real(Y)); axis('image'); title('a curvelet: spatial viewpoint');
subplot(1, 2, 2); colormap gray; imagesc(abs(F)); axis('image'); title('a curvelet: frequency viewpoint');

[SX, SY, FX, FY, NX, NY] = fdct_wrapping(C);
