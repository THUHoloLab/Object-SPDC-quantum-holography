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
%% Read hologram
image=im2double(rgb2gray(imread('./source/phase_object.jpg')));
back=im2double(rgb2gray(imread('./source/background.jpg')));

holo=image./back;
holo=holo(451:630,446:625);
figure(1),imshow(holo)
max(max(holo))

aver=mean(mean(holo));
holo=holo-aver+1;
figure(2),imshow(holo,[])

[m,n]=size(holo);
[M,N]=size(holo);


%% Experiment Parameters
wavelength = 810*10^(-9);% wavelength
pitch = 13*10^(-6);% pitch size of detector

for z=0.03:0.001:0.06
z
prop = Propagator_function(m, n, wavelength, pitch, z);% Transfer Function

%% Directly Reconstruction
measured=sqrt(holo);
recons = IFT((FT(holo)).*(prop));
pause(0.5)
figure(3),imshow(-abs(recons),[])
colormap(winter)
colorbar
clim([-1,-0.8])


%% Support
tic
support=zeros(m,n);
for ii=1:m
    for jj=1:n
        if (ii-(m/2)+5).^2+(jj-(n/2)-31).^2<23.^2
        support(ii,jj)=1;
        end
        if (ii-(m/2)+5).^2+(jj-(n/2)-79).^2<23.^2
        support(ii,jj)=1;
        end
        if (ii-(m/2)+5).^2+(jj-(n/2)+17).^2<23.^2
        support(ii,jj)=1;
        end
        if (ii-(m/2)+5).^2+(jj-(n/2)+65).^2<23.^2
        support(ii,jj)=1;
        end
    end
end

figure(4),imshow(support,[])



%% Iteration
Loops=100;
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

object=object.*support;
ph1=ph1.*support;


recons1=(1-object).*exp(1i.*ph1);

holo_field_updated = IFT((FT(recons1)).*conj(prop1));
A = abs(holo_field_updated);
phase = angle(holo_field_updated);
end

%%
    
A(start_m+1:start_m+m,start_n+1:start_n+n)=holo;
holo_field = A.*exp(1i.*phase);
recons1 = IFT((FT(holo_field)).*(prop1));
recons2=recons1(start_m+1:start_m+m,start_n+1:start_n+n);


%% Reconstruction

figure(5),imshow(abs(recons2),[]);
image_abs=abs(recons2);
image_abs=image_abs-min(min(image_abs));
image_abs=image_abs./max(max(image_abs));

figure(6),imshow(angle(recons2),[]);
image_phase=angle(recons2);
image_phase=image_phase-min(min(image_phase));
image_phase=image_phase./max(max(image_phase));


end



