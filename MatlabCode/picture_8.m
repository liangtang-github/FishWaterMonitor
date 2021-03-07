%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  药物浓度: 10TU,  1TU=2.3mg/L
%  开始时间：2017.7.31.9:08
%  加药时间：2017.7.31.13:09
%  结束时间：2017.8.1.8:36
%  按照一分钟一次：241分钟时下毒
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

dosingInv = 241;

h = figure;
set(h,'position',[300 150 900 580]);

% read 10TU TCP data
load('TCP10TU20170731.mat');
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
annotation(h,'textarrow',[0.187 0.213], [0.933 0.933],'String',{'13:09'});
xlabel('Sequentially');
ylabel('Signal intensity');
axis([0 size(IntegSignal,1) 0 1])
title('(a). 10TU 2,4,6-TCP exposure experiment original signal');

% EMD method
subplot(3,2,3);
set(gca,'position',[0.055,0.40,0.42,0.24]);
signal1 = emd(signal);
signal1 = (signal1 - min(signal1)) ./ (max(signal1) - min(signal1));
plot(signal1(:,4),'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.098 0.125], [0.288 0.288],'String',{'241'});
annotation(h,'textarrow',[0.164 0.132], [0.618 0.637],'String',{'258'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 250 500 750 1000 1250 1500]);
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
annotation(h,'textarrow',[0.580 0.607], [0.617 0.617],'String',{'241'});
annotation(h,'textarrow',[0.642 0.614], [0.618 0.638],'String',{'261'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 250 500 750 1000 1250 1500]);
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
annotation(h,'textarrow',[0.101 0.127], [0.619 0.619],'String',{'241'});
annotation(h,'textarrow',[0.163 0.136], [0.290 0.314],'String',{'276'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 250 500 750 1000 1250 1500]);
axis([0 size(signal3,1) 0 1])
title('(d). Wavelet method');

% unwinding method
subplot(3,2,6);
set(gca,'position',[0.53,0.0748,0.448,0.24]);
t=linspace(0,2*pi,length(signal));
L=15;
[F, err, S1, a, amp, ~, blaschke_z] = Unwinding_Blaschke(signal.',L,t);
[~, amplitude] = disp_uw(a,F);
wavesignal = amplitude(13,:);
wavesignal = (wavesignal - min(wavesignal)) ./ (max(wavesignal) - min(wavesignal));
signal4 = wavesignal;
plot(signal4,'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.580 0.607], [0.283 0.283],'String',{'241'});
annotation(h,'textarrow',[0.636 0.611], [0.293 0.313],'String',{'258'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 250 500 750 1000 1250 1500]);
axis([0 size(signal4,2) 0 1])
title('(e). Unwinding method');