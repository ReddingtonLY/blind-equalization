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