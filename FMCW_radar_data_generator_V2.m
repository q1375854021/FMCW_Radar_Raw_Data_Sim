c=3.0e8;             %Light Velocity.
B=1798.92e6;             %Bandwidth    change   Hz
K=29.982e12;             %Frequency Slop  change  Hz/S
T=B/K;                   %Ramp End Time  change S
Tc=160e-6;               %Chirp Cycle Time = Ramp End Time + Idle Time  =��100+60��e-6  change S
fs=10e6;                 %Sample rate change sps��How many times ADC sampling per second�� ��λsps  1sps=1Hz
f0=77e9;                 %Start Freq   change Hz
lambda=c/f0;             %Radar signal wavelength
d=lambda/2;              %Antenna array spacing   
nSamples = 512;      %ADC samples per chirp
allFrame = 3750;       %�ܹ��ɼ�����֡  ������汾���õ�allFrame��ʱ����mFFT
%%%
PRF = 250;%��λ�����Ƶ��
PRI=1/PRF;
Lsar = 14;
Nr = 512;%��ʱ������� 
Na = 3750;%��ʱ�������
V = Lsar/(PRI*Na);
%%%
fr = (0:Nr-1)/Nr*fs;
tr = (0:Nr-1)/fs;
Rr = c*fr/(2*K);
fa = (-Na/2:Na/2-1)/Na*PRF;
ta = (-Na/2:Na/2-1)*PRI;
Ra = V*ta;

%%%
%% ��Ŀ�����
target = [40 0 0];                  %Ŀ���λ��
radar = [0 0 0];
[target_number,~] = size(target);
sif = zeros(Na,Nr);          %��ƵƵ��

for j = 1:Na
    R = sqrt( (target(1)-radar(1)).^2 + (target(2)-(V*(ta(j)+tr)-radar(2))).^2);
    tau = 2*R/c;
    sif(j,:) = exp(1j*2*pi*K*tr.*tau + 1j*2*pi*f0*tau - 1j*pi*K*tau.^2);
end

%% ��Ŀ�����
% target = [10 0 0;
%           20 0 0;
%           30 0 0;
%           40 0 0];                  %Ŀ���λ��
% radar = [0 0 0];
% [target_number,~] = size(target);
% sif = zeros(Na,Nr);          %��ƵƵ��
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
xlabel('ʱ��');
title('��λ��ʱ��Ϊ0�ľ�������');

%���о��븵��Ҷ�任
fft_sif = fft(sif,Nr,2);
figure(2)
subplot(2,2,1)
mesh(fr,ta,abs(fft_sif));
xlabel('����Ƶ��');
ylabel('��λʱ��');
title('����fft��');
subplot(2,2,2)
plot(ta,abs(fft_sif(:,103)));
xlabel('��λʱ��');
title('������任��Ŀ�괦��λ���ΰ���');
subplot(2,2,3)
plot(ta,real(fft_sif(:,103)));
xlabel('��λʱ��');
title('������任��Ŀ�괦��λ����ʵ��');
subplot(2,2,4)
plot(ta,imag(fft_sif(:,103)));
xlabel('��λʱ��');
title('������任��Ŀ�괦��λ�����鲿');

H_RVP = exp(1j*pi*fr.^2/K);
fft_sif = fft_sif.*H_RVP;
figure(3)
subplot(2,2,1)
imagesc(fr,ta,abs(fft_sif));
xlabel('����Ƶ��');
ylabel('��λʱ��');
title('ȥ��RVP��');
subplot(2,2,2)
plot(ta,abs(fft_sif(:,103)));
xlabel('��λʱ��');
title('ȥ��RVP��Ŀ�괦��λ���ΰ���');
subplot(2,2,3)
plot(ta,real(fft_sif(:,103)));
xlabel('��λʱ��');
title('ȥ��RVP��Ŀ�괦��λ����ʵ��');
subplot(2,2,4)
plot(ta,imag(fft_sif(:,103)));
xlabel('��λʱ��');
title('ȥ��RVP��Ŀ�괦��λ�����鲿');


%���з�λ����Ҷ�任
fft_sif_2 = fftshift(fft(fft_sif,Na,1),1);
figure(4)
subplot(2,2,1)
mesh(fr,ta,abs(fft_sif_2));
xlabel('����Ƶ��');
ylabel('��λƵ��');
title('����fft��');
subplot(2,2,2)
plot(fa,abs(fft_sif_2(:,308)));
xlabel('��λƵ��');
title('��λ��任��Ŀ�괦��λ���ΰ���');
subplot(2,2,3)
plot(fa,real(fft_sif_2(:,308)));
xlabel('��λƵ��');
title('��λ��任��Ŀ�괦��λ����ʵ��');
subplot(2,2,4)
plot(fa,imag(fft_sif_2(:,308)));
xlabel('��λƵ��');
title('��λ��任��Ŀ�괦��λ�����鲿');

%% ƥ���˲�
for i=1:Nr                   
    Ka=2*V^2/(lambda*Rr(i));
    H = exp(1i*pi/Ka*fa.^2);
    fft_a_sif_RVP(:,i)= ifft(fft_sif_2(:,i).*H.',Na,1);
end
figure(5)
subplot(2,2,1)
imagesc(Rr,ta,abs(fft_a_sif_RVP));
xlabel('���������');
ylabel('��λʱ��');
title('��λ��ʱ��ѹ����');
subplot(2,2,2)
plot(ta,abs(fft_a_sif_RVP(:,308)));
xlabel('��λʱ��');
title('��λ��ʱ��ѹ����Ŀ�괦��λ���ΰ���');
subplot(2,2,3)
plot(ta,real(fft_a_sif_RVP(:,308)));
xlabel('��λʱ��');
title('��λ��ʱ��ѹ����Ŀ�괦��λ����ʵ��');
subplot(2,2,4)
plot(ta,imag(fft_a_sif_RVP(:,308)));
xlabel('��λʱ��');
title('��λ��ʱ��ѹ����Ŀ�괦��λ�����鲿');


