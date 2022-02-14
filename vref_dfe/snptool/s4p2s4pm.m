function s4pm = s4p2s4pm(s4p)
    M = [ 1  0 -1  0
          0  1  0 -1
	  1  0  1  0
	  0  1  0  1] / sqrt(2);
    invM = inv(M);
    nf   = size(s4p,1);
    s4pm = zeros(nf,4,4);
    for k = 1:nf
	s4pm(k,:,:) = shiftdim(M * shiftdim(s4p(k,:,:), 1) * invM, -1);
    end
end
