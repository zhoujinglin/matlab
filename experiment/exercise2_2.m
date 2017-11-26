clear; close all; clc;
H=fspecial('gaussian',9,2);
A=rgb2gray(imread('timg.jpg'));
AN=imnoise(A,'gaussian',0,0.02);
AN=double(AN);
 
TA=padarray(AN,[4,4],'replicate');
SFRes=zeros(size(AN));
for i=1:size(AN,1)
    for j=1:size(AN,2)
        for m=1:9
            for n=1:9
                SFRes(i,j)=SFRes(i,j)+TA(i+m-1,j+n-1)*H(m,n);
            end
        end
    end
end
%begin freq filt
[FH,f1,f2]=freqz2(H,size(AN));
FTR=fftshift(fft2(AN));
FTF=FTR.*FH;
FFRes=ifft2(ifftshift(FTF));
%end freq filt
FT=fftshift(fft2(double(A)));
FTSFR=fftshift(fft2(SFRes));
 
close all;
subplot(2,2,1),imshow(A);title('original');
subplot(2,2,2),imshow(mat2gray(abs(FT),[0,100000]));title('original FFT2');
subplot(2,2,3),imshow(mat2gray(AN));title('noisy');
subplot(2,2,4),imshow(mat2gray(abs(FTR),[0,100000]));title('noisy FFT2');
figure;
imshow(mat2gray(H,[0,0.05]));title('H');
figure;
subplot(2,3,1),imshow(mat2gray(SFRes));title('S Res');
subplot(2,3,2),imshow(mat2gray(abs(FTSFR),[0,100000]));title('S FFT2');
subplot(2,3,3),imshow(mat2gray(abs(FH)));title('FH Plane');
subplot(2,3,4),mesh(abs(FH));title('FH 3D');
subplot(2,3,5),imshow(mat2gray(FFRes));title('F Res');
subplot(2,3,6),imshow(mat2gray(abs(FTF),[0,100000]));title('F FFT2');
