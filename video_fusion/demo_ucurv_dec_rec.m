close all; clear all;
% 
disp(' ');
disp(' ');
disp(['In this demo, we test decomposition and reconstruction of 2-D ucurvelet']);
disp(['We can also estimate the running time']);
disp(' ');
r = input('Press <enter> to continue ...');

r = pi*[0.3 0.5 0.85 1.15];
alpha = 0.15;

N = [6 6 6 6];
lev = 2;

S = [512 512];
im = mkZonePlate(S);
% zeros data structure

% create curvelet window
F = ucurv_win(S, N, r, alpha);

tic
yg = ucurvdec(im, N, F);
toc

imr = ucurvrec(yg, N, F);

snr(im,imr)


