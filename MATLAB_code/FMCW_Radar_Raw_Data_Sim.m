%% 用于仿真FMCW雷达的中频信号
%相关的推导过程和注意事项，我会在相关的说明文档中指出。
clear;close all;clc;         
%% 配置SAR系统相关的参数
c=3.0e8;                 %光速   单位 m/s
B=1798.92e6;             %带宽   单位 Hz
K=29.982e12;             %调频率  单位 Hz/s
T=B/K;                   %调频时间 单位 s
% Tc=160e-6;               %Chirp Cycle Time = Ramp End Time + Idle Time=（100+60）e-6  change  S 应该暂时没有用  
fs=10e6;                 %距离向采样率 单位 Hz
f0=77e9;                 %起始频率 单位 Hz
lambda=c/f0;             %波长 单位 m
% d=lambda/2;              %天线的间距 单位 m  应该暂时无用
% nChirps = 12;             %一帧的chirp数，暂时无用
nSamples = 512;         %ADC在一个chirp内采了多少点
allFrame = 3750;        %一共采了多少帧，一帧设为一个chirp，也就是方位向点数
Lsar = 5;               %SAR平台运动长度
PRF = 250;              %方位向采样频率  也可认为是脉冲重复的频率
PRI=1/PRF;              %方位向采样间隔  就是chirp发射的间隔。


%% 计算相关的参数
Nr = nSamples;         %快时间，距离向采样点数
Na = allFrame;         %慢时间，方位向采样点数
V = Lsar/(PRI*(Na-1));     %SAR平台运动速度
fr = (0:Nr-1)/Nr*fs;   %距离向对应的频率
tr = (0:Nr-1)/fs;      %距离向对应的快时间
Rr = c*fr/(2*K);       %根据快时间计算距离向距离
fa = (-Na/2:Na/2-1)/Na*PRF;     %方位向对应的频率
ta = (-Na/2:Na/2-1)*PRI;        %方位向对应的慢时间
Ra = V*ta;                      %方位向对应的距离

%% 配置目标
% 格式[距离向距离，方位向距离，高度]
% 默认SAR平台的高度是0
target = [40 0 0;              %目标的位置
          41 1 0;
          44 -1 0;
          45 0 0;];                  
radar = [0 0 0];             %SAR平台的中心位置，实际上雷达应该略高一些
[target_number,~] = size(target);      %求取目标的个数
sif = zeros(Na,Nr);          %用于存储中频信号，列代表方位向，共有Na个数据。行代表距离向，共有Nr个数据
temp = zeros(Na,Nr);         %临时的数组，用于存储单个目标的回波

for i = 1:target_number      %对每个目标进行循环
    for j = 1:Na             %对每个方位向进行循环
        R = sqrt( (target(i,1)-radar(1)).^2 + (target(i,2)-(radar(2)+V*(ta(j)+tr))).^2 + (target(i,3)-radar(3)).^2 );     %求取目标和雷达之间的瞬时距离
        tau = 2*R/c;            %计算时延
        temp(j,:) = exp(1j*2*pi*K*tr.*tau + 1j*2*pi*f0*tau - 1j*pi*K*tau.^2);         %根据公式计算单目标回波
%         temp(j,:) = exp(-1j*2*pi*K*tr.*tau - 1j*2*pi*f0*tau + 1j*pi*K*tau.^2);      %两种均是变频结果，可根据实际情况使用         
    end
    sif = sif + temp;                  %对所有目标的回波进行叠加，即得到了最终的回波数据
end