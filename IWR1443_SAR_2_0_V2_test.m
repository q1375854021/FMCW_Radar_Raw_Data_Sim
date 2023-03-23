%%%This script is used to test the signal sampled by the IWR1443 ,if it's
%%%expressed by complex number.
%%%参考程序Kdil
%%%采样点大于420。
%仅使用距离向压缩和方位向压缩
clear all;close all;clc;
%% The parameter of IWR1443.
c=3.0e8;             %Light Velocity.
B=1798.92e6;             %Bandwidth    change   Hz
K=29.982e12;             %Frequency Slop  change  Hz/S
T=B/K;                   %Ramp End Time  change S
Tc=160e-6;               %Chirp Cycle Time = Ramp End Time + Idle Time  =（100+60）e-6  change S
fs=10e6;                 %Sample rate change sps（How many times ADC sampling per second） 单位sps  1sps=1Hz
f0=77e9;                 %Start Freq   change Hz
lambda=c/f0;             %Radar signal wavelength
d=lambda/2;              %Antenna array spacing   
nChirps = 12;
nSamples = 512;      %ADC samples per chirp
allFrame = 3750;       %总共采集多少帧  在这个版本中用的allFrame暂时代替mFFT
%%%
V = 0.09166667 ;           %单位 m/s 滑轨速度 5500mm/min  
Lsar = 1400e-3;
Nr = 512;%快时间采样点 
Na = 3750;%慢时间采样点
%%%

%%%
PRF = 250;%方位向采样频率
PRT=1/PRF;
% Ta=Na*PRT;
% Tr=Nr/fs;
% % fr=linspace(-fs/2,fs/2,Nr); 
% fr=linspace(0,fs,Nr); %%%去掉距离向FFT负频域
% ta=linspace(Ta,0,Na);  %慢时间   %感觉问题出在这个地方    %慢时间
% fa=linspace(-PRF/2,PRF/2,Na);
% tr=linspace(-Tr/2,Tr/2,Nr);%快时间 
fr = (0:Nr-1)./Nr.*fs;                  %距离向频率
tr = fr/K;                          %距离向时间
fa = (-Na/2:Na/2-1)./Na.*PRF;       %方位向多普勒频域
ta = (-Na/2:Na/2-1).*PRT;           %方位向时间
range_f_axis = fr;
range_r_axis = (c*range_f_axis)./(2*K);      %fft 快时间后距离坐标轴 
range_fang_axis = V.*ta;   %方位向坐标尺度换算 三种：滑轨运行的时间；天线长度；波束宽度        将慢时间变为中心对称，0~N-1变为 -N/2~N/2-1

d_res = c/(2*B);  %距离分辨率，参考PPT，用来确定x坐标轴

load sif;

%一行512点，对应快时间，一个Chirp的采样点
%一列3750点，对应慢时间，对应滑轨的位置
figure(1)
subplot(2,1,1)
plot(0:511,real(sif((fix(Na/2)),:)));  %Na为慢时间采样点数，也是3750  fix为朝0方向进行取整，Na/2，恰好为滑轨到中心的慢时间的快时间采样点
xlabel('快时间采样点');
ylabel('信号')
title('中频信号时域波形实数部分');
subplot(2,1,2)
plot(0:511,imag(sif((fix(Na/2)),:)));  %Na为慢时间采样点数，也是3750  fix为朝0方向进行取整，Na/2，恰好为滑轨到中心的慢时间的快时间采样点
xlabel('快时间采样点');
ylabel('信号')
title('中频信号时域波形虚数部分');


%% 首先对距离向进行FFT
nr = 0:Nr-1;        %range 距离向采样点
Rvp = exp(-1i*pi*K*(nr./fs).^2);
range_win = hamming(Nr);   %加个汉明窗   效果也不明显
for i = 1:Na
    sif(i,:) = (sif(i,:)-mean(sif(i,:)));   %剪掉直流，影响不大
%     sif(i,:) = fftshift(fft(fftshift(sif(i,:))));    %快时间fft
    sif(i,:) = fft(sif(i,:).*range_win').*Rvp;    %快时间fft
end

%确定距离坐标轴
% range_f_axis = ((0:Nr-1)-Nr/2)/Nr*fs;    %fft 快时间后频率坐标轴     fftshift永



figure(2)
imagesc(range_r_axis,range_fang_axis,abs(sif));
xlabel('距离向距离');
ylabel('方位向时间对应的距离');
title('距离向FFT后结果');

%% 方位向FFT
doppler_win = hamming(Na);        %减少截断引起的频谱泄露  
for ii=1:Nr    %快时间循环，对方位向进行处理  
             sif(:,ii) = (sif(:,ii)-mean(sif(:,ii)));       %信号维度3750*512   //一行数据代表一帧的chirp(1 Frame中第3个)采样点，共采了3750个
                                                            %一行512点，对应快时间，一个Chirp的采样点
                                                            %一列3750点，对应慢时间，对应滑轨的位置
                                                            %是想去掉均匀噪声吗，或者去掉直流分量
            sif(:,ii)=fftshift(fft((fftshift(sif(:,ii).*doppler_win))));    %在慢时间维度加窗，还进行了FFT
                                                            %一行512点，对应快时间，一个Chirp的采样点，
                                                            %一列3750点，对应慢时间的FFT，对应滑轨(方位向）FFT
end

%% 距离徙动矫正
Rd=c*fr/2/K;
delta_R = c.*fs./(2*K*Nr);
for i = 1:Na
    for j = 1:Nr
        delta_R_eta_feta = (lambda.^2.*Rd(j).*fa(i).^2)./(8*V.^2);
        delta_n = round(delta_R_eta_feta./delta_R);
        if(j+delta_n<=Nr)
            sif(i,j) = sif(i,j+delta_n);
        else
            sif(i,j) = 0;
        end
    end
end
figure(3)
imagesc(range_r_axis,range_fang_axis,abs(ifft(sif,Na)));
xlabel('距离向距离');
ylabel('方位向时间对应的距离');
title('距离徙动矫正后结果');


%% 检查sif
figure(4)
m = size(sif);
x = 1:m(2);
y = 1:m(1);
meshc(x,y,abs(sif));
xlabel('x  快时间  距离');
ylabel('y  慢时间  方位');

%% 匹配滤波后进行方位向IFFT
for i=1:Nr                   
%     Rd=-c*fr(ii)/2/K+Rc;                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Rd=-c*fr(i)/2/K;                
%     Kd=-2*V^2/lambda/Rd;
    Rd=c*fr(i)/2/K;                
    Kd=2*V^2/lambda/Rd;
    %sif(:,ii)=ifftshift(ifft(ifftshift(sif(:,ii).*exp(1i*pi*fa.^2/Kd).'.*exp(1i*4*pi/lamda*(Rd)))));
    sif(:,i)=ifftshift(ifft(ifftshift(sif(:,i).*exp(1i*pi*fa.^2/Kd).')));
end


% sif(abs(sif)>1700 & abs(sif)<6000)=sif(abs(sif)>1700 & abs(sif)<6000)*10;     %可以让目标点更明显，但是会掩盖掉这个范围外的目标
% sif(abs(sif)<1000 )=sif(abs(sif)<1000)*0;
% sif(:,220:511) = 0;
% sif(1:150,:) = 0;%%C

figure(5)
%% 检查sif
m = size(sif);
x = 1:m(2);
y = 1:m(1);
meshc(x,y,abs(sif));
xlabel('x  快时间  距离');
ylabel('y  慢时间  方位');


range_r_axis = (c*range_f_axis)./(2*K);      %fft 快时间后距离坐标轴 
range_a_axis = linspace(0,Lsar,Na) - Lsar/2;     %按照慢时间间隔进行分配

figure(6)
% imagesc(range_r_axis,range_a_axis,255-abs(sif));
imagesc(range_r_axis,range_a_axis,abs(sif));
axis equal
xlabel('距离向 m');
ylabel('方位向 m');
title('成像结果');
image_bad = sif;
save image_bad image_bad;
