function [inp,tnp,dnp,snp,f] = snp2impl(fin,snpin,varargin)

    prelen   = 0.040; % length of negative time against the entire period
    dcintlen = 0.1;   % length of integration to extrapolate DC gain if missing

    argk = 1;
    while (argk <= length(varargin))
	switch (varargin{argk})
	case 'prelen'
	    prelen = varargin{argk + 1};
	    argk = argk + 1;
	case 'dcintlen'
	    dcintlen = varargin{argk + 1};
	    argk = argk + 1;
	end
	argk = argk + 1;
    end

    % frequency delta
    df     = diff(fin(1:2));
    prepts = round(2 * fin(end) / df * prelen);	% number of data points in negative time

    % drop DC data point, if it has imaginary value
    if (fin(1) == 0)
	imDC = imag(snpin(1,:,:));
	if (sum(abs(imDC(:))) > 0)
	    fin   = fin(2:end);
	    snpin = snpin(2:end,:,:);
	end
    end

    % check if fin has DC
    if (abs(fin(1)) < 1e3)
	% DC is not missing
	f   = fin;
	snp = snpin;
    else
	% DC is missing

	% check if fin has low frequency entries
	if (round(fin(1) / df) == 1)
	    % LF entries are not missing

	    % extrapolate DC in time domain
	    f   = [0; fin];
	    dcg = snp2dcg(fin, snpin, prelen, dcintlen);
	    snp = [shiftdim(dcg, -1); snpin];

	else
	    % low frequency entries are missing

	    % extrapolate low frequency entries except DC in frequency domain
	    f = df * [0:round(fin(end)/df)]'; 
	    snpwithlf = interp1(fin, snpin, f(2:end), 'spline', 'extrap');

	    % extrapolate DC in time domain
	    dcg = snp2dcg(f(2:end), snpwithlf, prelen, dcintlen);
	    snp = [shiftdim(dcg, -1); snpwithlf];
	end
    end

    ts = 1 / 2 / f(end);	% time step

    [n0,n1,n2] = size(snp);
    n = n0 - 1;
    inp = zeros(2*n,n1,n2);
    tnp = zeros(2*n,n1,n2);
    dnp = zeros(n1,n2);
    for k1 = 1:n1
	for k2 = 1:n2
	    s1p  = snp(:,k1,k2);
	    phs  = angle(s1p);
	    dly  = (- phs(end) / pi) * ts;
	    s1p  = s1p .* exp(j * 2 * pi * f * dly);
	    s1psym = [s1p; conj(s1p(end-1:-1:2))];
	    i1psym = ifft(s1psym, 'symmetric');
%	    [pk, pkloc]  = max(abs(i1psym));
%	    sft  = mod(pkloc - prepts - 1, 2*n);
	    sft  = mod(1 - prepts - 1, 2*n);
%	    if (sft <= n + 1)
	    if (sft <= n)
		dly = dly + (sft - 0) * ts;
	    else
		dly = dly + (sft - 0 - 2*n) * ts;
	    end
	    inp(:,k1,k2) = i1psym([sft+1:end 1:sft]);
	    tnp(:,k1,k2) = [0:2*n-1]' * ts + dly;
	    dnp(k1,k2)   = dly;
	end
    end
end
