function s4p = s4pm2s4p(s4pm)
    M = [ 1  0 -1  0
          0  1  0 -1
	  1  0  1  0
	  0  1  0  1] / sqrt(2);
    invM = inv(M);
    nf   = size(s4pm,1);
    s4p  = zeros(nf,4,4);
    for k = 1:nf
	s4p(k,:,:) = shiftdim(invM * shiftdim(s4pm(k,:,:), 1) * M, -1);
    end
end
