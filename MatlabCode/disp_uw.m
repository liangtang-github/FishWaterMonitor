function [pha,amplitude] = disp_uw(a,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input : 'a' is a point in the unit disk;
%         'f' is real or analytic signal;
%
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = linspace(0,2*pi,length(f));
    phase_d = 0;
    amplitude = abs(f);
    for j = 1:length(a)
        phase_d = (1-abs(a(j))^2)./(1-2.*abs(a(j)).*cos(t-angle(a(j)))+ ...
            abs(a(j))^2)+phase_d;
    end
    Mphase_d = max(max(phase_d));
    phase_d = phase_d/Mphase_d;
    pha = phase_d*Mphase_d;
end