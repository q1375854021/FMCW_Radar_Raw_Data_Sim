clear all;close all;clc;
c=3.0e8;                 %Light Velocity.
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
%%���ںϳɿ׾����ȡ�
% �ϳɿ׾�����Խ��,Ŀ��ĺϳɿ׾�����Ҳ��Խ���������մ���Խ����λ�ֱ���Խ��
% �����Һϳɿ׾�����Ϊ1.4m�����߷�λ��׾��ֱȽ�խ���ᵼ���ӽ��ر��
% ֻ�н������Ŀ���п������������ĺϳɿ׾����̣���Զ����ĺϳɿ׾����Ⱥܿ��ܲ���������˽����Ŀ�귽λ��ֱ��ʻ�Ƚϸߡ�
% Զ��ķ�λ��ֱ��ʿ��ܻ�Ƚϵͣ�����ͨ����Lsar�ĳ�14��1.4��ֵ�������顣

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

target = [40 0 0;
          41 1 0;
          44 -1 0;
          45 0 0;];                  %Ŀ���λ��
radar = [0 0 0];
[target_number,~] = size(target);
sif = zeros(Na,Nr);          %��ƵƵ��
temp = zeros(Na,Nr);
% for i = 1:target_number
%     R0 = sqrt(target(i,1).^2 + target(i,2).^2);  
%     for m = 1:Na
%         R = sqrt(R0.^2 + (V*(ta(m)+tr)).^2);      %����
%         tau = 2*R/c;
%         temp(m,:) = exp(1j*2*pi*K*tr.*tau + 1j*2*pi*f0*tau - 1j*pi*K*tau.^2);
%     end
%     sif = sif + temp;
% end
%����һ�У�֮ǰ��Щ�����ţ������ٶ�û�б�ƽ�������·�λ��ѹ������ʧ�ܡ�
% �����д����ǰ��Ӧ���ڱ��Ƶ���ʽ����д���롣�����빫ʽ���Ӧ���Լ���д��ļ���

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
xlabel('ʱ��');
title('��λ��ʱ��Ϊ0�ľ�������');

%% ���о��븵��Ҷ�任
fft_sif = fft(sif,Nr,2);
figure()
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
% ���ھ�������Ҷ�任
% ������ѹ������ͨ����������Ҷ�任��ʵ�֣����ھ�������ʱ�����ѹ����Ӧ����һ��sinc�ͺ�����
% ����ֱ��ʿ����Լ��ơ�


%% ȥ��RVP
H_RVP = exp(1j*pi*fr.^2/K);
fft_sif = fft_sif.*H_RVP;
figure()
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
% ֮ǰ�о�RVP����ûʲôӰ�죬���ڿ���Ӱ�컹�ǽϴ�
% ���ԱȽ�RVP֮ǰ��֮��ľ���ѹ��Ч�������Կ���RVP�����Թ켣��������͵�һ��

% %% �����㶯����      %%�����в�ֵ��
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


%% �����㶯����      %��ֵ
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
%% �����㶯����      %��ֵ
for i = 1:Na
    for j = 1:Nr
        if(fr_new(i,j)<fr(1))      
            min(i) = j;          %��¼�����µ�������ǡ��С��ԭ���������ֵ
        elseif(fr_new(i,:)>fr(end))
            max(i) = j;          %��¼�����µ�������ǡ�ô���ԭ���������ֵ
            break;               %ΪʲôҪ��¼��Щ�أ���Ϊ��Щ�ᵼ���������ᳬ��ԭ��������ķ�Χ����ֵ��ǳ���׼ȷ��
        end                      %��˳�����Χ��������͸�ֵΪ0�Ϳ��ԡ� 
    end
end
for i = 1:Na
    sif_rcmc(i,:) = interp1(fr,fft_sif(i,:),fr_new(i,min(i):max(i)),'spline');
end

figure()
subplot(2,2,1)
imagesc(fr,ta,abs(sif_rcmc));
xlabel('����Ƶ��');
ylabel('��λƵ��');
title('�����㶯������');
subplot(2,2,2)
plot(fa,abs(sif_rcmc(:,308)));
xlabel('��λƵ��');
title('�����㶯����Ŀ�괦��λ���ΰ���');
subplot(2,2,3)
plot(fa,real(sif_rcmc(:,308)));
xlabel('��λƵ��');
title('�����㶯����Ŀ�괦��λ����ʵ��');
subplot(2,2,4)
plot(fa,imag(sif_rcmc(:,308)));
xlabel('��λƵ��');
title('�����㶯����Ŀ�괦��λ�����鲿');
% ��ͳ��RD�㷨�ǽ��з�λ����Ҷ�任����еģ������и�������һ���״��Ŀ��ľ�������̩�ս���չ���ɶ����
% ����Զ�࣬���Rrefѡ���������ļȿ��Դ���������������
% һ�����ǹ�ͨ�ġ�����������ʱ����һ�����������������㡣RrefҲ����ѡ�ڳ������ģ����Ƕ�ÿ�������Ŷ�����RCMC��
% ���Ͽ��Կ���������任��RD��Ҳû���˼������ƣ����ҹ�ʽ��ø����Ҳ�ֱ�ۣ�������ʱ��֮����о����㶯������
% ���ҽ��������������㣬����ֱ������ԭʼR��������㶯�����������f_rcmc�Ĺ�ʽ��ʾ��

%% ���з�λ����Ҷ�任
fft_sif_2 = fftshift(fft(sif_rcmc,Na,1),1);
figure()
subplot(2,2,1)
imagesc(fr,ta,abs(fft_sif_2));
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
% ��Ȼ���з�λ����Ҷ�任����Ҫ�ǹ۲�R-D���Ч��������Ҫ��Ŀ�Ļ��ǽ���ƥ���˲�����Ȼ���任ʱ����Ҳ�ǿ��Եġ�


%% ƥ���˲�
fft_a_sif_RVP = zeros(size(fft_sif));
for i=1:Nr                   
    Ka=2*V^2/(lambda*Rr(i));
    H = exp(1i*pi/Ka*fa.^2);
    fft_a_sif_RVP(:,i)= ifft(fft_sif_2(:,i).*H.',Na,1);
end
figure()
subplot(2,2,1)
mesh(Rr,ta,abs(fft_a_sif_RVP));
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


