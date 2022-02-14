function [wave,shiftcnt] = findStartIndex(wave, tpat)
    wave = wave(:);
    N    = length(wave);
    Npn  = length(tpat);
    Nosr = N / Npn;
    if Nosr == 1
        Nsub = 1;
    else
        Nsub = 4;
    end
    
    wavex2 = [wave; wave];
    acf1 = zeros(Npn*Nsub,1);
    for k = 1:Npn*Nsub
	Pmid = wavex2((k-1)*round(Nosr/Nsub)+[1:Nosr:N]);
	Pmid = Pmid - mean(Pmid);
%	acf1(k) = abs(mean(Pmid .* tpat));
	c = corrcoef([Pmid tpat.']);
	acf1(k) = abs(c(2,1));
    end

    [dmy,bestk1] = max(acf1);

    wvshift1 = circshift(wave, - (bestk1-1) * round(Nosr/Nsub));

    tpatX = filter(ones(Nosr,1),1,upsample(tpat,Nosr));

    acf2 = zeros(Nosr*2-1,1);
    for k = [1:Nosr*2-1]
	Pmid = circshift(wvshift1, k - Nosr);
	Pmid = Pmid - mean(Pmid);
%	acf2(k) = norm(Pmid .* tpatX);
	c = corrcoef([Pmid tpatX.']);
	acf2(k) = abs(c(2,1));
    end

    [dmy,bestk2] = max(acf2);

    wvshift2 = circshift(wvshift1, bestk2 - Nosr);

    wave = wvshift2;

    shiftcnt = - (bestk1-1) * round(Nosr/Nsub) + bestk2 - Nosr;
%% plot 
% plot(Pmid(1:Nosr*10));hold on
% grid on
% plot(tpatX(1:Nosr*10));


% function [ws,index] = findStartIndex(y,ideal_pulse)
% os_ui = 32;
% N = 8191;
% y_ = reshape(y,32,[]);vppy = mean(abs(y_(16,:)));
% vppi = mean(abs(ideal_pulse));
% scale = vppy/vppi;
% 
% ideal_pulse = ideal_pulse .* scale;
% % cross-correlation
% [coeff,lags] = xcorr(y(:),ideal_pulse);
% [~,id] = max(coeff);
% ws = id;% waveformshift right
% index = round(abs(lags(id))/os_ui)+1;
% % %debug
% ideal_pulseshift = circshift(ideal_pulse,(id));
% plot(ideal_pulseshift(32*1105:32*1125)*0.5);hold on
% grid on
% y_t = y(:);
% plot(y_t(32*1105:32*1125),'+');
% 
% %[~,exactIndexyasuo] = aligntpat(y,x); 
