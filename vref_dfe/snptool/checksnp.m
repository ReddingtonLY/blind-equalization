function [PQM, RQM, CQM, PWn, RWn, CMn, PMn, RMn, Rn] = checksnp(snp)

    Nf = size(snp, 1);
    N  = size(snp, 2);
    Ns = N * (N - 1);


    % check Passivity Quality Metric (PQM)

    A = 1.00001;
    B = 0.1;

    PMn = zeros(Nf, 1);
    for k = [1:Nf]
	Sk = shiftdim(snp(k,:,:), 1);
	PMn(k) = sqrt(max(eig(Sk'*Sk)));	% passivity measure at each frequency
    end
    PWn = max((PMn - A) / B, 0);		% passivity metric at each frequency
    PQM = max(100 / Nf * (Nf - sum(PWn)), 0);	% passivity quality metric

    % check Reciprocity Quality Metric (RQM)

    C = 0.000001;
    B = 0.01;

    RMn = zeros(Nf, 1);
    for k = [1:Nf]
	Sk = shiftdim(snp(k,:,:), 1);
	RMn(k) = sum(sum(abs(Sk - Sk'))) / Ns;	% reciprocity measure at each frequency
    end
    RWn = max((RMn - C) / B, 0);		% reciprocity metric at each frequency
    RQM = max(100 / Nf * (Nf - sum(RWn)), 0);	% reciprocity quality metric

    % check Causality Quality Metric (CQM)

    Vn = zeros(Nf-1, N*N);			% difference of SnP at two censecutive frequencies
    for k = [1:Nf-1]
	Vn(k,:) = reshape(snp(k+1,:,:) - snp(k,:,:), 1, N*N);
    end

    Rn  = zeros(Nf, 1);
    CMn = zeros(Nf, 1);
    for k = [1:Nf-2]
	Rn(k) = real(Vn(k+1)) .* imag(Vn(k)) - imag(Vn(k+1)) .* real(Vn(k));
	CMn(k) = Rn(k) / norm(Vn(k)) / norm(Vn(k+1));
    end
    CQM = max(100 * sum(max(Rn,0)) / abs(sum(Rn)), 0);

end

