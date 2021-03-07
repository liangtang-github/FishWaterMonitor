%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  药物名称：三氯酚,中华青鱂鱼
%  药物浓度:1TU,  1TU=2.3mg/L
%  开始时间：2017.4.19.10:15
%  加药时间：2017.4.19.14:15
%  结束时间：2017.4.20.14:15
%  按照一分钟一次：240分钟时下毒
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

dosingInv = 240;

h = figure;
set(h,'position',[300 150 900 580]);

% read 10TU TCP data
load('TCP1TU20170419.mat');
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
annotation(h,'textarrow',[0.161 0.187],[0.938 0.938],'String',{'14:15'});
xlabel('Sequentially');
ylabel('Signal intensity');
axis([0 size(IntegSignal,1) 0 1])
title('(a). 1TU 2,4,6-TCP exposure experiment original signal');

% EMD method
subplot(3,2,3);
set(gca,'position',[0.055,0.40,0.42,0.24]);
signal1 = emd(signal);
signal1 = (signal1 - min(signal1)) ./ (max(signal1) - min(signal1));
plot(signal1(:,4),'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.094 0.113],[0.619 0.617],'String',{'240'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
axis([0 size(signal1,1) 0 1])
xticks([0 250 500 750 1000 1250 1500]);
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
annotation(h,'textarrow',[0.571 0.592],[0.617 0.617],'String',{'240'});
annotation(h,'textarrow',[0.640 0.610],[0.616 0.638],'String',{'310'});
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
annotation(h,'textarrow',[0.094 0.112],[0.288 0.288],'String',{'240'});
annotation(h,'textarrow',[0.162 0.133],[0.294 0.315],'String',{'320'});
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
wavesignal = amplitude(14,:);
wavesignal = (wavesignal - min(wavesignal)) ./ (max(wavesignal) - min(wavesignal));
signal4 = wavesignal;
plot(signal4,'k','LineWidth',1.5);
hold on;
line([dosingInv,dosingInv],[0,1],'linestyle',':','LineWidth',2,'color','k');
grid on;
annotation(h,'textarrow',[0.571 0.592],[0.289 0.289],'String',{'240'});
annotation(h,'textarrow',[0.634 0.609],[0.293 0.313],'String',{'304'});
xlabel('Time span/(minute)');
ylabel('Signal intensity');
xticks([0 250 500 750 1000 1250 1500]);
axis([0 size(signal4,2) 0 1])
title('(e). Unwinding method');