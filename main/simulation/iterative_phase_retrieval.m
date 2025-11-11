% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Yiqian Yang, Liangcai Cao*
% %State Key Laboratory of Precision Measurement Technology and Instruments, Department of Precision Instruments, Tsinghua University, Beijing 100084, China
% %yang-yq22@mails.tsinghua.edu.cn
% %clc@tsinghua.edu.cn
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The code is written by Yiqian Yang 2025
% % The software is Matlab R2024a


clc;
clear all;
close all;
%% Generate Hologram

image=imread('./source/THU2.bmp');
image=im2double(image(:,:,1));
figure(1),imshow(image)


[m,n] = size(image);
[M,N] = size(image);


wavelength = 532*10^(-9); % wavelength
k = 2*pi/wavelength;
L0 = 5*10^(-3);
x0 = linspace(-L0/2,L0/2,m);
y0 = linspace(-L0/2,L0/2,m);
[x0,y0] = meshgrid(x0,y0);
z = 0.06;
F0 = exp(j*k*z)/(j*wavelength*z);
F1 = exp(j*k/2/z.*(x0.^2+y0.^2));
fF1 = fft2(F1);
fa1 = fft2(image);
Fuf1 = fa1.*fF1;
Uh = F0.*fftshift(ifft2(Fuf1));
Ih = Uh.*conj(Uh);    % image of hologram
figure(2),imshow(Ih,[0,max(max(Ih))/1])



%% Experiment Parameters
addpath('./Function');
pitch = 10*10^(-6);% pitch size of detector
prop = Propagator_function(m, n, wavelength, pitch, z);% Transfer Function

%% Directly Reconstruction
measured=sqrt(Ih);
recons = IFT((FT(Ih)).*(prop));
pause(0.5)
figure(3),imshow(abs(recons),[])
z


%% Support

support=zeros(m,n);

x = 200;  % x 坐标
y = 220;  % y 坐标
width = 90;  % 矩形宽度
height = 50;  % 矩形高度
support(y:y+height-1, x:x+width-1) = 1;



figure(4),imshow(support,[])

%% Iteration 
holo=Ih;
Loops=50;
M=1000;
N=1000;
A=ones(M,N);
phase=zeros(M,N);
start_m=(M-m)/2;
start_n=(N-n)/2;
prop1 = Propagator_function(M, N, wavelength, pitch, z);% Transfer Function
support=padarray(support,[(M-m)/2,(N-n)/2]);


for tt = 1:Loops
fprintf(': %d\n', tt)

A(start_m+1:start_m+m,start_n+1:start_n+n)=holo;
holo_field = A.*exp(1i.*phase);

recons1 = IFT((FT(holo_field)).*prop1);
object=1-abs(recons1);
ph1=angle(recons1);

for ii=1:M
    for jj=1:N
        if (object(ii,jj)<0)
            object(ii,jj)=0;
            ph1(ii,jj)=0;
        end
    end
end

object=object.*support;
ph1=ph1.*support;


recons1=(1-object).*exp(1i.*ph1);


% recons1 = 1-object;


holo_field_updated = IFT((FT(recons1)).*conj(prop1));
A = abs(holo_field_updated);
phase = angle(holo_field_updated);
end

%%
    
A(start_m+1:start_m+m,start_n+1:start_n+n)=holo;
holo_field = A.*exp(1i.*phase);
recons1 = IFT((FT(holo_field)).*(prop1));
recons2 = recons1(start_m+1:start_m+m,start_n+1:start_n+n);



%% Result
figure(5),imshow(abs(recons2),[]);
figure(6),imshow(angle(recons2),[]);
figure(7),imshow(abs(recons1),[]);
figure(8),imshow(angle(recons1),[]);