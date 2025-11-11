% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Yiqian Yang, Liangcai Cao*
% %State Key Laboratory of Precision Measurement Technology and Instruments, Department of Precision Instruments, Tsinghua University, Beijing 100084, China
% %yang-yq22@mails.tsinghua.edu.cn
% %clc@tsinghua.edu.cn
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The code is written by Yiqian Yang 2025
% % The software is Matlab R2024a


clc
clear all;
close all;
%% Generate Hologram

image=imread('./source/THU2.bmp');
image=im2double(image(:,:,1));
figure(1),imshow(image)


[m,n] = size(image);
[M,N] = size(image);


wavelength = 532*10^(-9);% wavelength
k = 2*pi/wavelength;


z = 0.02;
% Experiment Parameters
pitch = 5*10^(-6);% pitch size of detector
prop = Propagator_function(m, n, wavelength, pitch, z);% Transfer Function
Ih= IFT((FT(image)).*(prop));
Ih=abs(Ih).^2;
Ih=Ih./max(max(Ih));

figure(2),imshow(Ih,[0,max(max(Ih))/1])


variance = 0.1;   % noise
Ihn = imnoise(Ih, 'gaussian', 0, variance);
figure(3),imshow(abs(Ihn),[])


% Directly Reconstruction
measured=sqrt(Ihn);
recons = IFT((FT(Ihn)).*(prop));
figure(4),imshow(abs(recons),[])


PSNR = psnr(recons, image);
PSNR


function PSNR = psnr(f1, f2)
k = 8;
fmax = 2.^k - 1;
a = fmax.^2;
MSE =(double(im2uint8(f1)) -double( im2uint8(f2))).^2;
b = mean(mean(MSE));
PSNR = 10*log10(a/b);
end



