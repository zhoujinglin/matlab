clc,clear,close all;
I_A=[-3     0     3     3     0;
     -4     0     4     4     0;
     -4     0     4     4     0;
     -4     0     4     4     0;
     -3     0     3     3     0];
 [M, N]=size(I_A);
threshold=2; %设定阈值
I_A_threshold=zeros(M,N);
for i=1:M
    for j=1:N
        if I_A(i,j)>threshold
            I_A_threshold(i,j)=1;
        else
            I_A_threshold(i,j)=0;
        end
    end
end
subplot(1,2,1),imshow(I_A,[]),title('模板A');
subplot(1,2,2),imshow(I_A_threshold,[]),title('阈值分割图'); 