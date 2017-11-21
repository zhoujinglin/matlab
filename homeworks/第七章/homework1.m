clc,clear,close all;

I = [1, 1, 1, 0, 0; 1, 1, 1, 0, 0; 1, 1, 1, 0, 0; 1, 1, 1, 0, 0; 1, 1, 1, 0, 0];
A = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
B = [0, 1, 2; -1, 0, 1; -2, -1, 0];
I_A = imfilter(I, A, 'conv', 0, 'same');
I_B = imfilter(I, B, 'conv', 0, 'same');
subplot(1, 3, 1); imshow(I, []); title('I');
subplot(1, 3, 2); imshow(I_A, []); title('I与A卷积');
subplot(1, 3, 3); imshow(I_B, []); title('I与B卷积');
