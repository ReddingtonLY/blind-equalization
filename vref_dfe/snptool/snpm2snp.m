function snp = snpm2snp(snpm, posports, negports)
    N = size(snpm,2);
    M = eye(N);
    for k = [1:length(posports)]
	pos = posports(k);
	neg = negports(k);
	M(pos,neg) = -1;
	M(neg,pos) =  1;
    end
    M = M / sqrt(2);
    invM = inv(M);
    nf   = size(snp,1);
    snpm = zeros(nf,N,N);
    for k = 1:nf
	snpm(k,:,:) = shiftdim(invM * shiftdim(snp(k,:,:), 1) * M, -1);
    end
end
