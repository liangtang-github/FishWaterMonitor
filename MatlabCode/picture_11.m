%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  药物名称：铜离子（硫酸铜）,中华青鱂鱼
%  药物浓度: 1TU,  1TU=7.9mg/L
%  开始时间：2017.8.7.9:34
%  加药时间：2017.8.7.13:23
%  结束时间：2017.8.8.10:53
%  按照一分钟一次：229分钟时下毒
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

dosingInv = 229;

h = figure;
set(h,'position',[300 150 900 580]);

% read 10TU TCP data
load('CuSO41TU20170807.mat');
IntegSignal = VarName3 + VarName4 + VarName5 + VarName6 + VarName7 + ...
              VarName8 + VarName9 + VarName10;
IntegSignal = IntegSignal ./ 8;
signal = GetSignalByMinute(IntegSignal, 1200);

% orignal signal
subplot(3,2,[1, 2]);
set(gca,'position',[0.055,0.72,0.925,0.24]);
plot(IntegSignal,'k','LineWidth',1.5);
hold on;
line([dosingInv*1200,dosingInv*1200],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.154 0.126], [0.920 0.920],'String',{'13:23'});
xlabel('Sequentially');
ylabel('Signal intensity');
axis([0 size(IntegSignal,1) 0 1])
title('(a). 1TU CuSO_4 exposure experiment original signal');

% EMD method
subplot(3,2,3);
set(gca,'position',[0.055,0.40,0.42,0.24]);
signal1 = emd(signal);
signal1 = (signal1 - min(signal1)) ./ (max(signal1) - min(signal1));
plot(signal1(:,4),'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.113 0.086], [0.428 0.428],'String',{'229'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 600 1200 1800 2400 3000]);
axis([0 size(signal1,1) 0 1])
title('(b). EMD method');

% SG smooth method
subplot(3,2,4);
set(gca,'position',[0.53,0.40,0.448,0.24]);
signal2 = sgolayfilt(signal,3,61);
signal2 = (signal2 - min(signal2)) ./ (max(signal2) - min(signal2));
plot(signal2,'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.595 0.564], [0.465 0.464],'String',{'229'});
annotation(h,'textarrow',[0.789 0.815], [0.440 0.448],'String',{'1885'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 600 1200 1800 2400 3000]);
axis([0 size(signal2,1) 0 1])
title('(c). Savitzky Golay method');

% wavelet method
subplot(3,2,5);
set(gca,'position',[0.055,0.0748,0.42,0.24]);
[c,l] = wavedec(signal,3,'db3');
approx = appcoef(c,l,'db3');
approx = (approx - min(approx)) ./ (max(approx) - min(approx));
signal3 = interp(approx, 8);
plot(signal3,'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.115 0.087], [0.124 0.124],'String',{'229'});
annotation(h,'textarrow',[0.298 0.320], [0.107 0.132],'String',{'1895'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 600 1200 1800 2400 3000]);
axis([0 size(signal3,1) 0 1])
title('(d). Wavelet method');

% unwinding method
subplot(3,2,6);
set(gca,'position',[0.53,0.0748,0.448,0.24]);
t=linspace(0,2*pi,length(signal));
L=25;
[F, err, S1, a, amp, ~, blaschke_z] = Unwinding_Blaschke(signal.',L,t);
[~, amplitude] = disp_uw(a,F);
wavesignal = amplitude(17,:);
wavesignal = (wavesignal - min(wavesignal)) ./ (max(wavesignal) - min(wavesignal));
signal4 = wavesignal;
plot(signal4,'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.593 0.564], [0.193 0.191],'String',{'229'});
annotation(h,'textarrow',[0.890 0.865], [0.292 0.312],'String',{'2210'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 600 1200 1800 2400 3000]);
axis([0 size(signal4,2) 0 1])
title('(e). Unwinding method');