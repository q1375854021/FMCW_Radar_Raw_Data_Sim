c=3.0e8;             %Light Velocity.
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
%% 单目标测试
target = [40 0 0];                  %目标的位置
radar = [0 0 0];
[target_number,~] = size(target);
sif = zeros(Na,Nr);          %中频频率

for j = 1:Na
    R = sqrt( (target(1)-radar(1)).^2 + (target(2)-(V*(ta(j)+tr)-radar(2))).^2);
    tau = 2*R/c;
    sif(j,:) = exp(1j*2*pi*K*tr.*tau + 1j*2*pi*f0*tau - 1j*pi*K*tau.^2);
end

%% 多目标测试
% target = [10 0 0;
%           20 0 0;
%           30 0 0;
%           40 0 0];                  %目标的位置
% radar = [0 0 0];
% [target_number,~] = size(target);
% sif = zeros(Na,Nr);          %中频频率
% 
% for i = 1:target_number
%     R0 = sqrt(target(i,1).^2 + target(i,2).^2);
%     for j = 1:Na
%         R = sqrt(R0.^2 + V*(ta(j)+tr).^2);
%         tau = 2*R/c;
% %         sif(j,:) = sif(j,:) + exp(-1j*2*pi*f0*tau-1j*2*pi*K*tr.*tau + 1j*pi*K*tau.^2);
%         sif(j,:) = sif(j,:) + exp(1j*2*pi*K*tr.*tau + 1j*2*pi*f0*tau - 1j*pi*K*tau.^2);
%     end
% end



%%
save sif sif;
figure(1)
plot(tr,real(sif(Na/2,:)));
xlabel('时间');
title('方位向时间为0的距离向波形');

%进行距离傅里叶变换
fft_sif = fft(sif,Nr,2);
figure(2)
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

H_RVP = exp(1j*pi*fr.^2/K);
fft_sif = fft_sif.*H_RVP;
figure(3)
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


%进行方位向傅里叶变换
fft_sif_2 = fftshift(fft(fft_sif,Na,1),1);
figure(4)
subplot(2,2,1)
mesh(fr,ta,abs(fft_sif_2));
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

%% 匹配滤波
for i=1:Nr                   
    Ka=2*V^2/(lambda*Rr(i));
    H = exp(1i*pi/Ka*fa.^2);
    fft_a_sif_RVP(:,i)= ifft(fft_sif_2(:,i).*H.',Na,1);
end
figure(5)
subplot(2,2,1)
imagesc(Rr,ta,abs(fft_a_sif_RVP));
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


