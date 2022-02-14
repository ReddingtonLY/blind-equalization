function [best_ffecoeff,rocoeff,dfecoeff,best_decode,y_ui_s] = rxeq(param,algor)
% @Yuan 2021/4/1
% equalization to the waveform
Txfir = param.en_txfir;
% ideal pam signal
modE          = 4;                 % modulation level
bitpersym    = log2(modE);         % bits per symbols
tMode        = param.trainingMode;% if 0, using slicer on DFE, otherwise dk_hat = training seq
dk_hatMode   = param.dk_hatMode;  % if 0, use dk to update the tap weights, otherwise us dk_hat
if tMode == 0
    repeating = 1;
else 
    repeating = 1;	
end			  % 
repeat_times = 3; param.repeat_times = 3;
osr          = param.osr;
symRate      = param.bitRate/bitpersym;
nyfreq       = symRate/2;

table   = [-1 -1/3 1/3 1];	% pam 4
symbols = prbs13q(0)/3;		% prbs13q
%dk_hat  = symbols;

if Txfir == 1
    y_su = applyTxfir(param,symbols);  
else
    % without txfir
    for ii = 1:1:repeat_times
        y_su(:,(ii-1)*8191+1:8191*ii) = ones(osr,8191).*symbols;
    end
end
%% read s4p channel files
reads4p;
dt = mean(diff(tnp))

if param.pltLvl == 1
    figure;
    plot(fx*1e-9,db(abs(s21ddx)));hold on; grid on
    [~,id] = min(abs(nyfreq-fx));
    plot(fx(id)*1.e-9,20*log10(abs(s21ddx(id))),'o','LineWidth',4);
    sprintf('channel IL: %f dB',20*log10(abs(s21ddx(id))))   
end



[sbr,stp] = impl2sbr(tnp,inp,symRate*1e-9);
% ____________________________________________________________________________________ %
if  param.pltLvl == 1
    figure;plot(tnp,inp);
    grid on
    hold on;plot(tnp,sbr);
    xlabel('second');ylabel('voltage');%legend('impulse response','pulse response');hold off
end

% alignment 
len_sbr = floor(length(sbr)/osr); % ui
sbr = sbr(1:len_sbr*osr);
tpat = [zeros(15,1);1;zeros(len_sbr-16,1)]';
[sbr,shiftw] = findStartIndex(sbr,tpat);

osr          = param.osr;
if osr == 1
    %output = conv(y_su(16,:),inp(1:32:end),'same');
    output = filter(sbr(1:32:end),1,y_su(16,:));
    %output = conv(y_su(16,:),sbr(1:32:end),'same');
else
    output = conv(y_su(:),inp,'same');
    %output = filter(inp,1,y_su(:));
end

p2v = sbr(1:end-osr) - sbr(osr+1:end);
[~,idx0] = min(p2v);
[~,idx1] = max(p2v);
[~,idx2] = max(sbr);
zerophase  = interp1(p2v(idx0:idx1),[idx0:idx1],0,'spline');%cross zero
samphase = idx2 - zerophase;

N  = 8191*osr; % number of samples per repetition

[y_ui,shiftw] = findStartIndex(output(1:N),symbols);
%[y_ui,shiftw] = findStartIndex(output(end-N+1:end),symbols);

% add white nosie
y_ui = y_ui(:);
y_ui = reshape(awgn(y_ui,100,'measured'),osr,[]);

y_ui = reshape(y_ui,osr,[]);
% temp = y_ui;
% y_ui = [];
% y_ui = [zeros(32,10),temp,zeros(32,10)];
%% debug
% p = [15 1000 3  5  200 12 1       10       1];
% [y_avg,P,sigma_e,X,x] = get_pr(y_ui,symbols,p);
% figure;plot([0:1:length(P(:))-1]/32,P(:));xlabel('UI');ylabel('voltage');title('pulse response');

% setting
% ffe
ffeV = param.ffeRange;
c0   = param.mainCursor;
% dfe
Ndfe = param.Ndfe;% length of dfe
bmax = param.bmax;

[Nosr,patlen] = size(y_ui);best_SNRdb = -100;

if strcmpi(algor,'c-cma')||strcmpi(algor,'c-cma_t')
    
    muvec_ffe = [0.1];muvec_dfe = [0.1];
else
    muvec_ffe = [0.1];muvec_dfe = [0.01];
    
end


if param.iteTime ~= 0 
    for mui = 1:1:length(muvec_ffe)
        best_SNRdb = -100;
        if param.sampre == 0
            srange = 1:1:Nosr;
        else
            srange = param.sampre;
        end
        if repeating == 1
            subs = 0;
        else
            subs = ffeV;
        end
        dfe_snr = zeros(length(srange),patlen-subs);
        ffe_snr = zeros(length(srange),patlen-subs);
        yk_hat = zeros(1,patlen-subs);
        ynffeo = zeros(1,patlen-subs);
        dk = zeros(1,patlen-subs);
        dk_hat_s = cell(4,1);
        dk_hat_tgt = zeros(1,patlen-subs);
        if tMode == 1
            dk_hat = symbols;
        else
            dk_hat = zeros(1,patlen-subs);
        end
        error_dfe = zeros(1,patlen-subs);
        %error_ffe = zeros(1,patlen-subs);
        decoded_symbols = zeros(1,patlen-subs);
        for mid = srange
            % scan all possible phase
            % interl
            y_ui_s = y_ui;  
            %x_s = y_ui_s(mid,:);
            x_s = interp1([1:length(y_ui(:))],y_ui(:),[samphase:osr:patlen*32],'spline');
            [Nosr,patlen] = size(y_ui_s);
                  

        			% init
                    % ffe coefficient vector
                    ffeCoeff = zeros(1,ffeV); ffeCoeff(c0) = 1;
                    
                    clen = length(ffeCoeff);
                    % dfe coefficient vector
                    dfeCoeff = zeros(4,max(Ndfe,1));

                    % rls will update dfe as well
                    % cursor index
                    post = ffeV - c0;
                    pre = c0 - 1;


                    %----------------------------------------------------------------------------------------%
                    %---------------------------internal parameters set--------------------------------------%
                    %----------------------------------------------------------------------------------------%
                    % variable step size (not applied)
                    q = 2;
                    p = 2;

                    mu_ffe = muvec_ffe(mui);          % step size
                    mu_dfe = muvec_dfe(mui);


                    %----------------------------------------------------------------------------------------%
                    %----------------------------------------------------------------------------------------%
                    % [1 2 3 4 5 6 7 8] --> [6 7 8 1 2 3 4 5 6 7 8] xn_e
                    % based on repeating pattern property
                    if repeating == 1
                        [xn_e] = expX(x_s,post,pre,patlen,clen);
                    else
                        xn_e = x_s;
                    end
                    % target
                    vpp = mean(abs(xn_e))*1.5;%;max(sbr)
                    %----------------------------------------------------------------------------------------%
                    %----------------------------------------------------------------------------------------%

                    %% start of iteration
                    % errf stand for error function for each algorithm

                    for ite = 1:param.iteTime

                        for n = 1:1:patlen-subs
%                             if n+(ite-1)*(patlen-subs) < 20000
%                                 Ndfe = 0;
%                             else
%                                 Ndfe = 1;
%                             end
                            %-------------------------------------------------------------------------------------%
        					% X1 is ffe input(inversed)
                            X1 = xn_e(n+clen-1:-1:n)';% xn_e(n+clen-1:-1:n)
                            if ~iscolumn(X1)
                                X1 = X1';
                            end
                            % ynffeo is ffe output
                            ynffeo(n) = ffeCoeff*X1;
                            %-------------------------------------------------------------------------------------%
                            %------------------------------------------DFE----------------------------------------%
                            %-------------------------------------------------------------------------------------%
                           
                            vref_dfe_process
                            
                            % dfe error and ffe error
                            if tMode == 0
                                alphbet   = table.*vpp;
                                [~,idx]   = min(abs(dk(n)-alphbet));
                                dk_hat(n) = table(idx);decoded_symbols(n) = dk_hat(n);
                            elseif tMode == 1
                                decoded_symbols(n) = dk_hat(n);
                            end
                            
                         
                            error_dfe(n) = dk(n) - dk_hat(n)*vpp;
                            
                            %vpp_ffeo = mean(abs(ynffeo))/2;
                            %error_ffe(n) = dk(n) - dk_hat(n)*vpp_ffeo;

                            % for snr
                            %error_dfe_1(n) = dk(n) - dk_hat(n)*vpp_dfe;


                            %----------------------------------------------------------------------------------------%
                            %------------------------------------ffe-coefficient-update------------------------------%
                            %----------------------------------------------------------------------------------------%
                            if strcmpi(algor,'c-cma_T')
                                R  =  vpp;%mean(abs(dk))*1.5;%max(sbr);
                                R1 = R;R2 = R/3;

                                a1 = (R1+R2)/2; a2 = (R1-R2)/2;

        						errf(n)  = a2-abs(a1-abs(dk(n)));

                                M(n)     = sign(a1-abs(dk(n)))*sign(dk(n));
                                %                         R  =  mean(abs(ynffeo))*1.5;%max(sbr);
                                %                         R1 = R;R2 = R/3;
                                %
                                %                         a1 = (R1+R2)/2; a2 = (R1-R2)/2;
                                %
                                % 						errf(n)  = a2-abs(a1-abs(ynffeo(n)));
                                %
                                %                         M(n)     = sign(a1-abs(ynffeo(n)))*sign(ynffeo(n));

                                ffeCoeff = ffeCoeff - (2 * mu_ffe * errf(n) * M(n) * X1)' ;

                                %**************************************************************************************%
                                %**************************************************************************************%
                                %                         R  = vpp;%max(sbr); %mean(abs(dk))*1.5;
                                %                         R1 = R;R2 = R/3;
                                %
                                %                         a1 = (R1+R2)/2; a2 = (R1-R2)/2;
                                %
                                % 						errf(n)  = a2-abs(a1-abs(dk(n)));
                                %
                                %                         M(n)     = sign(a1-abs(dk(n)))*sign(dk(n));
                                % dfe
                                for ii=1:1:min(Ndfe,n-1)

                                    if 0%dk_hatMode == 1
                                        temp_tar = dk_hat(n-ii)*vpp;
                                    else
                                        temp_tar = dk(n-ii);
                                    end
                                    if  n == 1
                                        update_c = 4;
                                    else
                                        switch(dk_hat(n-1))
                                            case 1
                                                update_c = 1;
                                            case 1/3
                                                update_c = 2;
                                            case -1/3
                                                update_c = 3;
                                            case -1
                                                update_c = 4;
                                        end
                                    end

                                    if (abs(dfeCoeff(update_c,ii))< R*bmax)
                                        dfeCoeff(update_c,ii) = dfeCoeff(update_c,ii) + 2 * mu_dfe * errf(n) * M(n) * temp_tar;
                                        %dfeCoeff(ii) = dfeCoeff(ii) + mu_dfe * error_dfe(n) * temp_tar;
                                    else
                                        dfeCoeff(update_c,ii) = sign(dfeCoeff(update_c,ii)) * R * bmax;
                                    end
                                end

                                %----------------------------------------------------------------------------------------%
                                %----------------------------------------------------------------------------------------%
                            elseif strcmpi(algor,'c-cma')
                                % p = 2, q = 1
                                R  = max(sbr); %mean(abs(dk))*1.5;%vpp;
                                R1 = R;R2 = R/3;

                                a1 = (R1^2+R2^2)/2; a2 = (R1^2-R2^2)/2;
                                theta = sqrt(a1+a2);

        						errf(n)  = dk(n)^2 - a1 - a2;

                                ffeCoeff = ffeCoeff - (mu_ffe * errf(n) * (X1*(dk(n)-theta)*(dk(n)+theta)^2 + X1*(dk(n)+theta)*(dk(n)-theta)^2))' ;

                                % dfe
                                for ii=1:1:min(Ndfe,n-1)

                                    if dk_hatMode == 1
                                        temp_tar = dk_hat(n-ii)*R;
                                    else
                                        temp_tar = dk(n-ii);
                                    end
                                    if  n == 1
                                        update_c = 4;
                                    else
                                        switch(dk_hat(n-1))
                                            case 1
                                                update_c = 1;
                                            case 1/3
                                                update_c = 2;
                                            case -1/3
                                                update_c = 3;
                                            case -1
                                                update_c = 4;
                                        end
                                    end

                                    if (abs(dfeCoeff(update_c,ii))< R*bmax)
                                        %dfeCoeff(ii) = dfeCoeff(ii) + mu_dfe * errf(n) * (dk(n-ii)*(dk(n)-theta)*(dk(n)+theta)^2 + dk(n-ii)*(dk(n)+theta)*(dk(n)-theta)^2) ;
                                        dfeCoeff(update_c,ii) = dfeCoeff(update_c,ii) + mu_dfe * error_dfe(n) * temp_tar;
                                    else
                                        dfeCoeff(update_c,ii) = sign(dfeCoeff(update_c,ii)) * R * bmax;
                                    end
                                end

                                %----------------------------------------------------------------------------------------%
                                %----------------------------------------------------------------------------------------%
                                % LMS %
                            elseif strcmpi(algor,'lms')

                                ffeCoeff = ffeCoeff - (2 * mu_ffe*(error_dfe(n))* X1)';

                                if dk_hatMode == 1
                                    temp_tar = dk_hat*vpp;
                                else
                                    temp_tar = dk;
                                end
                                if n == 1
                                    % assuming previous symbol -1
                                    dfeCoeff(4,:) = dfeupdate(dfeCoeff(4,:),Ndfe,n,vpp,bmax,mu_dfe,temp_tar,error_dfe);
                                else
                                    switch(dk_hat(n-1))
                                        case 1
                                            dfeCoeff(1,:) = dfeupdate(dfeCoeff(1,:),Ndfe,n,vpp,bmax,mu_dfe,temp_tar,error_dfe);
                                        case 1/3
                                            dfeCoeff(2,:) = dfeupdate(dfeCoeff(2,:),Ndfe,n,vpp,bmax,mu_dfe,temp_tar,error_dfe);
                                        case -1/3
                                            dfeCoeff(3,:) = dfeupdate(dfeCoeff(3,:),Ndfe,n,vpp,bmax,mu_dfe,temp_tar,error_dfe);
                                        case -1
                                            dfeCoeff(4,:) = dfeupdate(dfeCoeff(4,:),Ndfe,n,vpp,bmax,mu_dfe,temp_tar,error_dfe);
                                    end
                                end
                                    
                                %----------------------------------------------------------------------------------------%
                                %----------------------------------------------------------------------------------------%
                            else
                                error('no algorithm found!');
                            end
                            %----------------------------------------------------------------------------------------%
                            %----------------------------------------------------------------------------------------%
                            % calcaulate snr after * iteration
                            tt = n+(ite-1)*(patlen-subs);
                            if  tt > 1
                                SNRdb_ffe = calSNR(ynffeo,param);
%                                 for ii = 1:1:4
%                                     SNRdb_dfes(ii) = calSNR(dk_hat_s{1},param);
%                                 end
%                                 SNRdb_dfe = sum(SNRdb_dfes)/4;        
                                SNRdb_dfe = calSNR(dk,param);
                                dfe_snr(mid,n+(ite-1)*(patlen-ffeV)) = SNRdb_dfe;
                                ffe_snr(mid,n+(ite-1)*(patlen-ffeV)) = SNRdb_ffe;
                                %SNRdb_dfe = db(vpp_dfe/(norm(error_dfe_1)/sqrt(length(error_dfe_1))))-2;
                                %dfe_snr(mid+1,n+(ite-1)*(patlen-clen)) = SNRdb_dfe;
                                %ffe_snr(ins+1,n+(ite-1)*(patlen-clen)) = SNRdb;
                                %avg_eh(n+(ite-1)*(patlen-clen)) = eyeh(dk);
                            end
%                             tap0(n+(ite-1)*(patlen-clen)) = vpp;
%                             vffeo(n+(ite-1)*(patlen-clen)) = vpp_ffeo;
%                             dfe_m(:,n+(ite-1)*(patlen-clen)) = dfeCoeff;
%                             ffe_m(:,n+(ite-1)*(patlen-clen)) = ffeCoeff;
%                             err_mo_dfe(n+(ite-1)*(patlen-clen)) = error_dfe(n);
%                             err_mo_ffe(n+(ite-1)*(patlen-clen)) = error_ffe(n);
                        end
                        if ite == 1
                            errf = error_dfe;
                            error1ffe = errf;
                            error1dfe = error_dfe;
                        end
                        % AGC
                        %dk = dk * vpp / (mean(abs(dk)) * 1.5);

                    end % iter loop

                    if SNRdb_dfe > best_SNRdb
                        best_sam = mid;
                        best_SNRdb = SNRdb_dfe;
                        best_ffecoeff = ffeCoeff;
                        best_dfecoeff = dfeCoeff;
                        best_yn = ynffeo;
                        best_dk = dk;
                        best_dfe_snr = dfe_snr(mid,:);
                        best_ffe_snr = ffe_snr(mid,:);
                        best_decode = decoded_symbols; % normalized decoded symbols
                    end                

        end % mid loop

        % show iteration vs snr
        %best_ffe_snr_new(mui,:) = best_ffe_snr(201:end);
        best_dfe_snr_new(mui,:) = best_dfe_snr(201:end);
        % reset for different step size
        if length(muvec_ffe) > 1
            if mui == length(muvec_ffe)
                best_dfe_snr = [];
                best_ffe_snr = [];
                ynffeo = [];
                dk = [];
            end
        end

    end % mu loop


    % output
    rocoeff = round(best_ffecoeff/sum(abs(best_ffecoeff))*63)
    best_dfecoeff'
    vpp
    dfecoeff = best_dfecoeff./vpp



    if param.pltLvl == 1
    %     figure;
    %     plot(avg_eh(201:end));title('average eye height');grid on
    %     
    %     figure;
    %     plot(best_dk,'*');title('slicer input');hold on
    %     plot(dk_hat.*vpp_dfe,'*');hold off
    %     
    %     figure
    %     plot(best_yn,'*');title('ffe output');hold on
    %     plot(yk_hat.*vpp_ffeo,'*');hold off
    %     
	% dfe SNR with different step size
    % 	figure;
    %     colormap hot;
    %     cmap = colormap;
    %     for mui = 1:1:length(muvec_ffe)
    %         Plot_color = cmap(mui*20,:);  
    %         plot(best_dfe_snr_new(mui,:),'linewidth',2,'Color',Plot_color);hold on
    %         legstr{mui} = sprintf('mu = %2d',muvec_ffe(mui));
    %         
    %     end
    % 
    %     grid on
    % 	legend(legstr);xlabel('iteration times');ylabel('SNR dB');
	
	% compare ffe SNR and DFE SNR
	figure; 
	plot(best_dfe_snr_new);
	%plot(best_ffe_snr_new(1,:));grid on
	%legend('dfe','ffe');
    xlabel('iteration times');ylabel('SNR dB');
	
	% TAP weight changing
    % 	figure;
    % 	for ii = 1:length(dfeCoeff)
    % 		plot(dfe_m(ii,:));hold on
    % 	end
    % 	legend('tap1','tap2','tap3');title('dfe tap weight changes')
    % 	grid on;hold off
    %     
    %     % ffe tap weight changes
    % 	figure;
    % 	for ii = 1:length(ffeCoeff)
    % 		plot(ffe_m(ii,:));hold on
    % 	end
    % 	title('ffe tap weight changes')
    %     
    %     % vpp changes
    %     figure;
    % 	plot(tap0);hold on
    %     plot(vffeo);
    %     legend('dk','yffeo');
    %     title('vpp changes');grid on
    %     
    %     % error
    %     figure;
    %     plot(err_mo_dfe);title('dfe error');
    %     figure;
    %     plot(err_mo_ffe);title('ffe error');
    end

    % debug
    best_c = 0;
    for ii = 1:1:8191
    c = corrcoef(circshift(yk_hat,ii-1),symbols); 
    if c(2,1)>best_c
        best_c = c(2,1);
        best_s = ii-1;
    end
    end    
    
else
    % pr mode
    y_ui_s = y_ui;
    best_decode = dk_hat;
    rocoeff = 0;
    dfecoeff = 0;
    best_ffecoeff = 0;
end

end


%-----------------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------------%
function [wave,shiftcnt] = aligntpat(wave, tpat)
wave = wave(:);
N    = length(wave);
Npn  = length(tpat);
Nosr = N / Npn;

Nsub = 4;

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
end

%-----------------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------------%
function [avgEh] = eyeh(x)
% the minium eye height of the half signal  


l3 = mean(abs(x))*1.5; l2 = l3/3; l1 = -l2; l0 = -l3;
sl1 = (l3+l2)/2; sl2 = (l0+l1)/2;
xl3 = x(find(x > sl1));
xl2 = x(find(x < sl1 & x > 0));
xl1 = x(find(x > sl2 & x < 0));
xl0 = x(find(x < sl2));

N0 = round(length(xl0)/2); % half 
N1 = round(length(xl1)/2); % half 
N2 = round(length(xl2)/2); % half 
N3 = round(length(xl3)/2); % half 

if N0 ~= 0 && N1 ~= 0 && N2 ~=0 && N3 ~= 0
    eye3 = min(xl3(end-N3:end))-max(xl2(end-N2:end));
    eye2 = min(xl2(end-N2:end))-max(xl1(end-N1:end));
    eye1 = min(xl1(end-N1:end))-max(xl0(end-N0:end));
else
    eye1 = 0;
    eye2 = 0;
    eye3 = 0;
end

avgEh = (eye3+eye2+eye1)/3;

end
%-----------------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------------%
function [SNRdb] = calSNR(ynffeo,param)

if strcmpi(param.mode,'pam4')
    l3 = mean(abs(ynffeo))*1.5; l2 = l3/3; l1 = -l2; l0 = -l3;
    sl1 = (l3+l2)/2; sl2 = (l0+l1)/2;
    % SNR after ffe equalization
    snr3 = mean(ynffeo(find(ynffeo > sl1)).^2)/var(ynffeo(find(ynffeo > sl1)));
    snr2 = mean(ynffeo(find(ynffeo < sl1 & ynffeo > 0)).^2)/var(ynffeo(find(ynffeo < sl1 & ynffeo > 0)));
    snr1 = mean(ynffeo(find(ynffeo > sl2 & ynffeo < 0)).^2)/var(ynffeo(find(ynffeo > sl2 & ynffeo < 0)));
    snr0 = mean(ynffeo(find(ynffeo < sl2)).^2)/var(ynffeo(find(ynffeo < sl2)));
    s = [snr0,snr1,snr2,snr3];
    s(isnan(s)) = [];
    %snr = min([snr0,snr1,snr2,snr3]);
    snr = mean(s);
    SNRdb = 10*log10(snr);
    
    if isnan(snr3)||isnan(snr2)||isnan(snr1)||isnan(snr0)
        SNRdb = -inf;
    end
elseif strcmpi(param.mode,'nrz')
    
    % SNR after ffe equalization
    snr1 = mean(ynffeo(find(ynffeo > 0)).^2)/var(ynffeo(find(ynffeo >  0)));
    snr0 = mean(ynffeo(find(ynffeo < 0)).^2)/var(ynffeo(find(ynffeo < 0)));
    s = [snr0,snr1];
    s(isnan(s)) = [];
    %snr = min([snr0,snr1,snr2,snr3]);
    snr = mean(s);
    SNRdb = 10*log10(snr);
    
    if isnan(snr1)||isnan(snr0)
        SNRdb = -inf;
    end
end
end
%-----------------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------------%
function [dfeCoeff] = dfeupdate(dfeCoeff,Ndfe,n,vpp,bmax,mu_dfe,dk_hat,err)
% dfe coefficient update
for ii=1:1:min(Ndfe,n-1)
    if (abs(dfeCoeff(ii))< vpp*bmax)
        dfeCoeff(ii) = dfeCoeff(ii) + (2*mu_dfe*err(n)*dk_hat(n-ii)')';
    else
        dfeCoeff(ii) = sign(dfeCoeff(ii))*vpp*bmax;
    end
end
end
%-----------------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------------%
function [norm_tx_coeff,Rmt]=findtxcoeff(m,param,P_t,rk)   
     % IEEE802.3ck 162.9.3.1.1
uilen = param.osr;
Np = param.Np;

cm = zeros(length(m),param.numPre+param.numPost+1);  % 162-2
em = zeros(length(m),1);  % 162-4
Rm = zeros(uilen*Np,param.numPre+param.numPost+1,length(m));
for m_idx = 1:length(m)
    
    for jj = 1:1:uilen*Np
        for ii = -param.numPre:1:param.numPost
            mm = m(m_idx);
            if (mm+jj-uilen*ii)>=1 && (mm+jj-uilen*ii)<=uilen*Np
                Rm(jj,ii+param.numPre+1,m_idx) = rk(mm+jj-uilen*ii);
            else
                Rm(jj,ii+param.numPre+1,m_idx) = 0;
            end
        end
    end
    
    cm(m_idx,:) = ((Rm(:,:,m_idx)'*Rm(:,:,m_idx))^-1)*Rm(:,:,m_idx)'*P_t(1:uilen*Np);
    em(m_idx) = sum((P_t(1:uilen*Np)-(Rm(:,:,m_idx)*cm(m_idx,:)')).^2);
end
[~,min_em] = min(em);
norm_tx_coeff = cm(min_em,:);

Rmt = Rm(:,:,min_em);
end
%-----------------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------------%

function [xn_e] = expX(x_s,post,pre,patlen,clen)

xn_e = zeros(1,patlen+clen-1);

if post ~= 0
    xn_e(1:post) = x_s(end-post+1:end);
end

xn_e(post+1:end-pre) = x_s;

if pre ~= 0
    xn_e(end-pre+1:end) = x_s(1:pre);
end



end