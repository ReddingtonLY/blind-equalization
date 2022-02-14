function s2pm = s2p2s2pm(s2p)
    M = [ 1  -1
	  1   1] / sqrt(2);
    invM = inv(M);
    nf   = size(s2p,1);
    s2pm = zeros(nf,2,2);
    for k = 1:nf
	s2pm(k,:,:) = shiftdim(M * shiftdim(s2p(k,:,:), 1) * invM, -1);
    end
end
