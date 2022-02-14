function s2p = s2pm2s2p(s2pm)
    M = [ 1  -1
	  1   1] / sqrt(2);
    invM = inv(M);
    nf   = size(s2pm,1);
    s2p  = zeros(nf,2,2);
    for k = 1:nf
	s2p(k,:,:) = shiftdim(invM * shiftdim(s2pm(k,:,:), 1) * M, -1);
    end
end
