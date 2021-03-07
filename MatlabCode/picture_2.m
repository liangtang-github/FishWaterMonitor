clear;
clc;
close all;

h = figure;
set(h,'position',[300 200 900 400]);

% Behavioral signals in a safe state.
load('TCP10TU20170731.mat');
IntegSignal = VarName3 + VarName4 + VarName5 + VarName6 + VarName7 + ...
              VarName8 + VarName9 + VarName10;
IntegSignal = IntegSignal ./ 8;
saftySignal = GetSignalByMinute(IntegSignal, 1200);
saftySignal = saftySignal(1:200);
saftySignal = (saftySignal - min(saftySignal)) ./ (max(saftySignal) - min(saftySignal));
subplot(2,2,1);
set(gca,'position',[0.055,0.59,0.42,0.35])
plot(saftySignal,'k','LineWidth',1.5);
grid on;
yticks([0 0.25 0.50 0.75 1.0]);
xticks([0 25 50 75 100 125 150 175 200]);
axis([0 size(saftySignal,1) 0 1])
xlabel('Time span/(minute)');
ylabel('Behavior intensity');
title('(a). Safty state behavior signal');


% Biological clock behavior signal.
load('Nontoxic20170419.mat');
IntegSignal = VarName3 + VarName4 + VarName5 + VarName6 + VarName7 + ...
              VarName8 + VarName9 + VarName10;
IntegSignal = IntegSignal ./ 8;
bioSignal = GetSignalByMinute(IntegSignal, 1200);
bioSignal = (bioSignal - min(bioSignal)) ./ (max(bioSignal) - min(bioSignal));
subplot(2,2,2);
set(gca,'position',[0.53,0.59,0.456,0.35])
plot(bioSignal,'k','LineWidth',1.5);
grid on;
yticks([0 0.25 0.50 0.75 1.0]);
xticks([0 500 1000 1500 2000 2500 3000 3500]);
axis([0 size(bioSignal,1) 0 1])
xlabel('Time span/(minute)');
ylabel('Behavior intensity');
title('(b). Biological clock behavior signal');


% Abnormal signals in TCP
load('TCP10TU20170731.mat');
IntegSignal = VarName3 + VarName4 + VarName5 + VarName6 + VarName7 + ...
              VarName8 + VarName9 + VarName10;
IntegSignal = IntegSignal ./ 8;
tcpSignal = GetSignalByMinute(IntegSignal, 1200);
tcpSignal = (tcpSignal - min(tcpSignal)) ./ (max(tcpSignal) - min(tcpSignal));
subplot(2,2,3);
set(gca,'position',[0.055,0.095,0.42,0.35])
plot(tcpSignal,'k','LineWidth',1.5);
grid on;
yticks([0 0.25 0.50 0.75 1.0]);
xticks([0 200 400 600 800 1000 1200 1400]);
axis([0 size(tcpSignal,1) 0 1])
xlabel('Time span/(minute)');
ylabel('Behavior intensity');
title('(c). Abnormal behavior signal in TCP');


% Abnormal signals in CuSO4
load('CuSO410TU20170803.mat');
IntegSignal = VarName3 + VarName4 + VarName5 + VarName6 + VarName7 + ...
              VarName8 + VarName9 + VarName10;
IntegSignal = IntegSignal ./ 8;
cusoSignal = GetSignalByMinute(IntegSignal, 1200);
cusoSignal = (cusoSignal - min(cusoSignal)) ./ (max(cusoSignal) - min(cusoSignal));
subplot(2,2,4);
set(gca,'position',[0.53,0.095,0.456,0.35])
plot(cusoSignal,'k','LineWidth',1.5);
grid on;
yticks([0 0.25 0.50 0.75 1.0]);
xticks([0 50 100 150 200 250 300 350 400]);
axis([0 size(cusoSignal,1) 0 1])
xlabel('Time span/(minute)');
ylabel('Behavior intensity');
title('(d). Abnormal behavior signal in CuSO4');







