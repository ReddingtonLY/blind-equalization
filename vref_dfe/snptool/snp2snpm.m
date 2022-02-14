function snpm = snp2snpm(snp, posports, negports)
    N = size(snp,2);
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
	snpm(k,:,:) = shiftdim(M * shiftdim(snp(k,:,:), 1) * invM, -1);
    end
end
