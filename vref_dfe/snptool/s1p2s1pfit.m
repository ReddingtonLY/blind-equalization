function [s1pfit, delay] = s1p2s1pfit(fin, s1p, varargin)

    % fit S1P frequency response

    fout       = fin;
    nfin       = length(fin);
    nfout      = nfin;
    nfit       = nfin;
    maxlossdB  = 20;
    weight     = ones(nfin,1);
    unitDCgain = false;

    argk = 1;
    while (argk <= length(varargin))
	switch(varargin{argk})
	case 'fout'
	    fout  = varargin{argk + 1};
	    nfout = length(fout);
	    argk  = argk + 1;
	case 'maxnfit'
	    nfit = varargin{argk + 1};
	    argk = argk + 1;
	case 'maxlossdB'
	    maxlossdB = varargin{argk + 1};
	    argk = argk + 1;
	case 'weight'
	    weight = varargin{argk + 1};
	    argk = argk + 1;
	case 'unitDCgain'
	    unitDCgain = true;
	end
	argk = argk + 1;
    end

    % check number of data points to fit

    nfitmax = find(db(s1p(1:nfit)) < db(s1p(1)) - maxlossdB, 1);

    if (length(nfitmax) ~= 0)
	nfit = nfitmax;
    end

    % check delay in time domain

    [i1p, t1p, d1p] = snp2impl(fin, s1p);

    [pk, pkidx] = max(i1p);
    idxofs = find(i1p(pkidx:-1:1) < 0, 1);	% offset to last negative value before peak
    if (length(idxofs) > 0)
	negidx = pkidx - idxofs + 1;		% index of last negative value before peak
	posidx = negidx + 1;			% index of first positive value after negidx
	i1pneg = i1p(negidx);
	i1ppos = i1p(posidx);
	t1pneg = t1p(negidx);
	t1ppos = t1p(posidx);
	delay  = (t1pneg * i1ppos - t1ppos * i1pneg) / (i1ppos - i1pneg);	% zero crossing time
    else
	delay  = 0;	% no delay
    end

    % fit in frequency domain

    if (unitDCgain)
	Xin     = [sqrt(fin)  fin];
	Xout    = [sqrt(fout) fout];
	Xweight = [weight     weight];
    else
	Xin     = [ones(nfin,1)  sqrt(fin)  fin];
	Xout    = [ones(nfout,1) sqrt(fout) fout];
	Xweight = [weight        weight     weight];
    end

    coeff  = (Xin(1:nfit,:) .* Xweight(1:nfit,:)) \ (db(s1p(1:nfit)) .* weight(1:nfit));

    % generate output S1P as fitting result

    gain   = 10 .^ (Xout * coeff / 20);
    phase  = snp2minph(fout, gain);
    s1pfit = gain .* exp(j * phase) .* sec2snp(fout, delay);

end
