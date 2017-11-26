clc,clear,close all;

x = importdata('seismic_251_301_2ms.txt');
[p,q] = size(x);
fs = 500;
r = 1;
dt = 1/fs;
I = (p-1)*dt;
t0 = 1.8;
t0_r = t0+I*(r-1);
s1 = x(:,r);
n = 0:p-1;
t = t0_r+n/fs;
%振幅随时间变化曲线
figure(1);plot(t,s1);
xlabel('t/ms');ylabel('振幅');
title('第一列振幅随时间变化曲线');
%频谱中心化
y = fft(s1,p-1);
mag = abs(y);
y = fftshift(y);
mag0 = abs(y);
M = length(y);
%中心频率图/有效频谱
f = (0:M-1)*fs/M;
f0 = f-f(M/2);
figure(2);
subplot(2,1,1),plot(f0,mag0);
xlabel('频率/Hz');ylabel('振幅');
title('全部频率');grid on;
subplot(2,1,2),plot(f(1:M/2),mag(1:M/2));
xlabel('频率/Hz');ylabel('振幅');
title('有效频率');grid on;
%生成Butterworth矩阵
H = zeros(3,M);
JS = [1 2 4];
D0 = 100;
for n = 1:M
  H(1,n) = 1/(1+((n-M/2)/D0)^(2*JS(1)));
  H(2,n) = 1/(1+((n-M/2)/D0)^(2*JS(2)));
  H(3,n) = 1/(1+((n-M/2)/D0)^(2*JS(3)));
end
%画出butterworth函数
figure(3);
subplot(3,1,1),plot(f0,abs(H(1,:)),'-b');hold on;
subplot(3,1,2),plot(f0,abs(H(2,:)),'-b');hold on;
subplot(3,1,3),plot(f0,abs(H(3,:)),'-b');hold on;
xlim([f0(1),f0(M)]);xlabel('f0/Hz');

z1 = y'.*H(2,:);
figure(4),plot(f0,abs(z1)),grid on;
