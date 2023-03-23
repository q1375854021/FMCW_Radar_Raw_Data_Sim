clear all;close all;clc;
c=3.0e8;                 %Light Velocity.
B=1798.92e6;             %Bandwidth    change   Hz
K=29.982e12;             %Frequency Slop  change  Hz/S
T=B/K;                   %Ramp End Time  change S
Tc=160e-6;               %Chirp Cycle Time = Ramp End Time + Idle Time  =（100+60）e-6  change S
fs=10e6;                 %Sample rate change sps（How many times ADC sampling per second） 单位sps  1sps=1Hz
f0=77e9;                 %Start Freq   change Hz
lambda=c/f0;             %Radar signal wavelength
d=lambda/2;              %Antenna array spacing   
nSamples = 512;      %ADC samples per chirp
allFrame = 3750;       %总共采集多少帧  在这个版本中用的allFrame暂时代替mFFT
%%%
PRF = 250;%方位向采样频率
PRI=1/PRF;
Lsar = 14;         
%%关于合成孔径长度。
% 合成孔径长度越长,目标的合成孔径历程也就越长，多普勒带宽越宽，方位分辨率越宽。
% 假设我合成孔径长度为1.4m，天线方位向孔径又比较窄，会导致视角特别宽。
% 只有近距离的目标有可能走完完整的合成孔径历程，而远距离的合成孔径长度很可能不完整，因此近距的目标方位向分辨率会比较高。
% 远距的方位向分辨率可能会比较低，可以通过将Lsar改成14，1.4等值进行试验。

Nr = 512;%快时间采样点 
Na = 3750;%慢时间采样点
V = Lsar/(PRI*Na);
%%%
fr = (0:Nr-1)/Nr*fs;
tr = (0:Nr-1)/fs;
Rr = c*fr/(2*K);
fa = (-Na/2:Na/2-1)/Na*PRF;
ta = (-Na/2:Na/2-1)*PRI;
Ra = V*ta;

%%%

target = [40 0 0;
          41 1 0;
          44 -1 0;
          45 0 0;];                  %目标的位置
radar = [0 0 0];
[target_number,~] = size(target);
sif = zeros(Na,Nr);          %中频频率
temp = zeros(Na,Nr);
% for i = 1:target_number
%     R0 = sqrt(target(i,1).^2 + target(i,2).^2);  
%     for m = 1:Na
%         R = sqrt(R0.^2 + (V*(ta(m)+tr)).^2);      %本行
%         tau = 2*R/c;
%         temp(m,:) = exp(1j*2*pi*K*tr.*tau + 1j*2*pi*f0*tau - 1j*pi*K*tau.^2);
%     end
%     sif = sif + temp;
% end
%在这一行，之前少些了括号，导致速度没有被平方，导致方位向压缩总是失败。
% 因此在写代码前，应当在边推导公式，编写代码。代码与公式相对应可以减少写错的几率

for i = 1:target_number
    for j = 1:Na
        R = sqrt( (target(i,1)-radar(1)).^2 + (target(i,2)-(V*(ta(j)+tr)-radar(2))).^2);
        tau = 2*R/c;
        temp(j,:) = exp(1j*2*pi*K*tr.*tau + 1j*2*pi*f0*tau - 1j*pi*K*tau.^2);
    end
    sif = sif + temp;
end


save sif sif;
figure()
plot(tr,real(sif(Na/2,:)));
xlabel('时间');
title('方位向时间为0的距离向波形');

%% 进行距离傅里叶变换
fft_sif = fft(sif,Nr,2);
figure()
subplot(2,2,1)
mesh(fr,ta,abs(fft_sif));
xlabel('距离频率');
ylabel('方位时间');
title('距离fft后');
subplot(2,2,2)
plot(ta,abs(fft_sif(:,103)));
xlabel('方位时间');
title('距离向变换后目标处方位向波形包络');
subplot(2,2,3)
plot(ta,real(fft_sif(:,103)));
xlabel('方位时间');
title('距离向变换后目标处方位向波形实部');
subplot(2,2,4)
plot(ta,imag(fft_sif(:,103)));
xlabel('方位时间');
title('距离向变换后目标处方位向波形虚部');
% 关于距离向傅里叶变换
% 距离向压缩可以通过距离向傅里叶变换来实现，由于距离向有时宽，因此压缩后应当是一个sinc型函数。
% 距离分辨率可以自己推。


%% 去除RVP
H_RVP = exp(1j*pi*fr.^2/K);
fft_sif = fft_sif.*H_RVP;
figure()
subplot(2,2,1)
imagesc(fr,ta,abs(fft_sif));
xlabel('距离频率');
ylabel('方位时间');
title('去除RVP后');
subplot(2,2,2)
plot(ta,abs(fft_sif(:,103)));
xlabel('方位时间');
title('去除RVP后目标处方位向波形包络');
subplot(2,2,3)
plot(ta,real(fft_sif(:,103)));
xlabel('方位时间');
title('去除RVP后目标处方位向波形实部');
subplot(2,2,4)
plot(ta,imag(fft_sif(:,103)));
xlabel('方位时间');
title('去除RVP后目标处方位向波形虚部');
% 之前感觉RVP可能没什么影响，现在看来影响还是较大。
% 可以比较RVP之前和之后的距离压缩效果，可以看见RVP后明显轨迹变得鲜明和单一。

% %% 距离徙动矫正      %%不进行插值，
% delta_f = fs/Nr;
% delta_show = zeros(Na,Nr);
% for i = 1:Na
%     for j = 2:Nr
%         f_rcmc = 2*K/c*(sqrt(Rr(j).^2 + (V*ta(i)).^2)-Rr(j));
%         delta = round(f_rcmc/delta_f);
%         delta_show(i,j) = delta;
%         sif_rcmc(i,j) = fft_sif(i,mod(j+delta,Nr)+1);
%     end
% end


%% 距离徙动矫正      %插值
delta_f = fs/Nr;
sif_rcmc = zeros(size(fft_sif));
f_rcmc = zeros(size(fft_sif));
for i = 1:Na
    for j = 2:Nr
        f_rcmc(i,j) = 2*K/c*(sqrt(Rr(j).^2 + (V*ta(i)).^2)-Rr(j));
    end
end
fr_new = ones(length(fa),1)*fr+f_rcmc;
min = ones(size(fr_new));
max = ones(size(min)).*length(fr);
%% 距离徙动矫正      %插值
for i = 1:Na
    for j = 1:Nr
        if(fr_new(i,j)<fr(1))      
            min(i) = j;          %记录下来新的坐标轴恰好小于原来坐标轴的值
        elseif(fr_new(i,:)>fr(end))
            max(i) = j;          %记录下来新的坐标轴恰好大于原来坐标轴的值
            break;               %为什么要记录这些呢？因为这些会导致新坐标轴超出原有坐标轴的范围，插值会非常不准确。
        end                      %因此超出范围的坐标轴就赋值为0就可以。 
    end
end
for i = 1:Na
    sif_rcmc(i,:) = interp1(fr,fft_sif(i,:),fr_new(i,min(i):max(i)),'spline');
end

figure()
subplot(2,2,1)
imagesc(fr,ta,abs(sif_rcmc));
xlabel('距离频率');
ylabel('方位频率');
title('距离徙动矫正后');
subplot(2,2,2)
plot(fa,abs(sif_rcmc(:,308)));
xlabel('方位频率');
title('距离徙动矫正目标处方位向波形包络');
subplot(2,2,3)
plot(fa,real(sif_rcmc(:,308)));
xlabel('方位频率');
title('距离徙动矫正目标处方位向波形实部');
subplot(2,2,4)
plot(fa,imag(sif_rcmc(:,308)));
xlabel('方位频率');
title('距离徙动矫正目标处方位向波形虚部');
% 传统的RD算法是进行方位向傅里叶变换后进行的，但那有个条件，一是雷达和目标的距离利用泰勒近似展开成二次项。
% 二是远距，因此Rref选到成像中心既可以代替整个成像区域。
% 一二条是共通的。比如近距成像时，第一条近似条件不再满足。Rref也不能选在成像中心，而是对每个距离门都进行RCMC。
% 综上可以看到，将其变换到RD域也没有了计算优势，而且公式变得复杂且不直观，不如在时域之间进行距离徙动矫正。
% 而且近似条件不再满足，不如直接利用原始R计算距离徙动矫正，如计算f_rcmc的公式所示。

%% 进行方位向傅里叶变换
fft_sif_2 = fftshift(fft(sif_rcmc,Na,1),1);
figure()
subplot(2,2,1)
imagesc(fr,ta,abs(fft_sif_2));
xlabel('距离频率');
ylabel('方位频率');
title('距离fft后');
subplot(2,2,2)
plot(fa,abs(fft_sif_2(:,308)));
xlabel('方位频率');
title('方位向变换后目标处方位向波形包络');
subplot(2,2,3)
plot(fa,real(fft_sif_2(:,308)));
xlabel('方位频率');
title('方位向变换后目标处方位向波形实部');
subplot(2,2,4)
plot(fa,imag(fft_sif_2(:,308)));
xlabel('方位频率');
title('方位向变换后目标处方位向波形虚部');
% 仍然进行方位向傅里叶变换，主要是观察R-D域的效果，最主要的目的还是进行匹配滤波，当然不变换时域卷积也是可以的。


%% 匹配滤波
fft_a_sif_RVP = zeros(size(fft_sif));
for i=1:Nr                   
    Ka=2*V^2/(lambda*Rr(i));
    H = exp(1i*pi/Ka*fa.^2);
    fft_a_sif_RVP(:,i)= ifft(fft_sif_2(:,i).*H.',Na,1);
end
figure()
subplot(2,2,1)
mesh(Rr,ta,abs(fft_a_sif_RVP));
xlabel('距离向距离');
ylabel('方位时间');
title('方位向时域压缩后');
subplot(2,2,2)
plot(ta,abs(fft_a_sif_RVP(:,308)));
xlabel('方位时间');
title('方位向时域压缩后目标处方位向波形包络');
subplot(2,2,3)
plot(ta,real(fft_a_sif_RVP(:,308)));
xlabel('方位时间');
title('方位向时域压缩后目标处方位向波形实部');
subplot(2,2,4)
plot(ta,imag(fft_a_sif_RVP(:,308)));
xlabel('方位时间');
title('方位向时域压缩后目标处方位向波形虚部');


