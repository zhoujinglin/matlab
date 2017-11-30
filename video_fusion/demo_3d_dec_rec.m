close all; clear all;
% 
disp(' ');
disp(' ');
disp(['In this demo, we demonstrate that 3-D ucurvelet']);
disp(['is perfect reconstruction']);
disp(' ');

rr = input('Press <enter> to continue ...');

tic

close all; clear all;
% size of data cube
Sz = [64 32 128];
% configuration , the number of direction at each size of pyramid
% total #direction is 3*Cf^2

Cf = [3 6];

% uncomment to have different number of direction for different dimension.
% Cf = [6 3 12;3 3 12;6 6 24];

% default window parameter.
r = pi*[0.3 0.5 0.85 1.15];
alpha = 0.15;

% F2 = ucurvwin3d(Sz*[1 1 1], Cf, r, alpha);

F2 = ucurvwin3d(Sz, Cf, r, alpha);

im = rand(Sz);

ydec = ucurvdec3d(im, Cf, F2);

imr = ucurvrec3d(ydec, Cf, F2);

snr(im,imr)
% ---------------------------------------------------------
% 
disp(' ');
disp(' ');
disp(['Now we the 3-D ucurvelet sparse window version']);
disp(['This version can run for large data cube']);
disp(' ');

rr = input('Press <enter> to continue ...');

Sz = 384*[1 1 1];
% configuration , the number of direction at each size of pyramid
% total #direction is 3*Cf^2

% Cf = [3 3];

Cf = [6 12;6 6;6 6];

% default window parameter.
r = pi*[0.3 0.5 0.85 1.15];
alpha = 0.15;

im = single(rand(Sz));

tic
% 3 sec for size 64 
% 161 sec for size 256
[F2, ind, cf] = ucurvwin3d_s(Sz, Cf, r, alpha);
toc

tic
% 1 sec for size 64
% 98 sec for size 256
ydec = ucurvdec3d_s(im, Cf, F2, ind, cf );
toc
% pack;
tic
imr = ucurvrec3d_s(ydec, Cf, F2, ind, cf );
toc
max(abs(im(:)-imr(:)))

snr(im,imr)


% ---------------------------------------------------------

tic
ydec = ucurvdec3d(im, Cf, r);
toc

tic
imr = ucurvrec3d(ydec, Cf, r);
toc

snr(im,imr)

