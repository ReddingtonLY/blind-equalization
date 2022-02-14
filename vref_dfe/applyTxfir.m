% apply Txfir to ideal waveform
% always set osr = 32 in this function
function [y_su] = applyTxfir(param,symbols) 
    osr          = 32;  repeat_times = param.repeat_times;
    % 802.3ck 162.9.3.1.3 Coefficient initialization preset3
    % c(-3) c(-2)  c(-1)  c(0)  c(1)
    %   0     0   -0.075  0.75   0
    % 162.9.3.1.1 Linear fit to the measured waveform Rm(j,i+4) MNp by 5
    % matrix
    if mod(osr,2) == 0
        m = -osr/2:1:osr/2-1;
    else
        m = -(osr-1)/2:1:(osr-1)/2-1;
    end
    param.Np = 200; param.Dp = 15;
    param.numPre  = 3;%param.mainCursor-1;
    param.numPost = 1;%param.ffeRange - param.mainCursor;
    
    P_I = zeros(param.Np*osr,1);P_I(15*osr+1:16*osr) = ones(osr,1);
    [norm_tx_coeff,Rmt] = findtxcoeff(m,param,P_I,P_I);
    %% temp test-------------------------------%
    txfir = param.txfir;%[-0.0104    0.0532   -0.1824    0.4645   -0.2894];
    txfir = txfir/sum(abs(txfir));
    %------------------------------------------%
    PR = reshape(txfir*Rmt',osr,[]);
    PR(:,end+1) = zeros(osr,1);
    %% linear fit use normalized x(n) as in 3ck 162.9.3.1.1
    xr = circshift(symbols,-param.Dp);  %left shift
    X = zeros(param.Np,length(symbols));
    X(1,:) = xr;
    %% X1
    for ii = 1:1:param.Np-1
        X(ii+1,:) = circshift(xr,ii);
    end
    X(end+1,:)=ones(1,length(symbols));
    out = PR*X;
    for ii = 1:1:repeat_times
        y_su(:,(ii-1)*8191+1:8191*ii) = ones(osr,8191).*out;
    end
  

end

%-----------------------------------------------------------------------------------------------------------%
%-----------------------------------------------------------------------------------------------------------%
function [norm_tx_coeff,Rmt]=findtxcoeff(m,param,P_t,rk)   
     % IEEE802.3ck 162.9.3.1.1
uilen = 32;
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