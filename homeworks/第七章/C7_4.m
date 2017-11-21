clc,clear,close all;
image=imread('C7_4.jpg');
image_th=(image>80);
imshow(image_th);