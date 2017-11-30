close all; clear all;
% 

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
disp(' ');
disp(['Estimate basis of UCURV 2-D version A']);
disp(['We can also estimate the norm of 2-D curvelet']);
disp(' ');
r = input('Press <enter> to continue ...');
% size or resulting image
S = 256;

% default UDCT parameters
alpha = 0.15;
r = pi*[0.3 0.5 0.85 1.15];

% direction configuration, changed to see different #direction or
% resolution
N = [6 12 24];

% create curvelet window
F = ucurv_win(S, N, r, alpha);

im256 = zeros(S);
% zeros data structure
yg = ucurvdec(im256, N, F);

lev = 2;
nlev = N(lev-1);

% insert one into subband at appropriate place
% mostly horizontal band
in1 = 1;
cs = size(yg{lev}{in1*nlev/2});
st1 = cs(1)/2;
st2 = cs(2)/(nlev/2);
for in2 = 1:(nlev/2);
    yg{lev}{(in1-1)*(nlev/2)+in2}(st1/2+(in1-1)*st1, round(0.5*st2+(in2-1)*st2)) = 0+1*j;
end


% mostly vertical band
in1 = 2;
cs = size(yg{lev}{in1*nlev/2});

if nlev == 6
    st1 = cs(1)/2;
    st2 = cs(2)/(nlev/2);

    for in2 = 1:(nlev/2);
        yg{lev}{(in1-1)*(nlev/2)+in2}(st1/2+(in1-1)*st1, round(0.5*st2+(in2-1)*st2)) = 0+1*j;
    end
else
    st1 = 0.5*cs(1)/(nlev/4);
    st2 = cs(2)/2;

    for inr = 1:nlev/4
        for inc = 1:2
            yg{lev}{(nlev/2)+(inc-1)*nlev/4+inr}(round(0.5*cs(1)+st1/2+(inr-1)*st1), round(0.5*st2+(inc-1)*st2)) = 0+1*j;
        end
    end
end


imr = ucurvrec(yg, N, F);

% create white border ------------------------------------------------

cs = size(imr);
st1 = cs(1)/2;
st2 = cs(2)/(nlev/2);

wht = 0.8* max(imr(:));
% middle white
imr(st1,:) = wht;

imr(1:st1,round(st2:st2:(S-1))) = wht;

% second half
if nlev == 6
    imr(st1+1:(S-1),round(st2:st2:(S-1))) = wht;
else
    imr(st1+1:(S-1),cs(2)/2) = wht;
    imr(round(st1+1:st1/(nlev/4):(S-1)),:) = wht;
end   
figure;    
imagesc(imr);axis off; axis equal


% --------------------------------------------------------------------





