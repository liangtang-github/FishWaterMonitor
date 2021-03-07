function [dealSignal] = GetSignalByMinuteOriginal(signal, thrd)
    len = size(signal, 1);
    dealSignal = zeros(floor(len ./ thrd), 1);
    
    for i = 1:size(dealSignal,1)
        dealSignal(i,1) = sum(signal((i-1)*thrd+1:i*thrd, 1)) ./ thrd;
    end
end