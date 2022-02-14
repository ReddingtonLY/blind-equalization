function [snpout,dnp,fint,snpint] = extrapsnp(fin,snpin,fout,varargin)

    % interpolate S parameter in different frequency

    method = 'makima';		% method of interpolation to step response in time domain
    force  = 'forcep';		% enforce given S parameter gain and phase in freq domain
    zerof  = 'nozerof';		% enforce zero for >x2 of given frequency
    causal = 'nocausal';	% enforce causality in time-domain
    maxge  = 10;		% max gain enforcement in dB
    extlim = 1.4;		% extrapolation limit of each step
    suckdb = 20;		% depth of suck out in dB

    k = 1;
    while (k <= length(varargin))
	opt = varargin{k};
	switch(opt)
	case {'spline','linear','makima'}
	    method = opt;
	case {'force','forceg','noforce'}
	    force  = opt;
	case {'zerof','nozerof'}
	    zerof = opt;
	case {'causal','nocausal'}
	    causal = opt;
	case 'maxforceg'
	    k = k + 1;
	    maxge = varargin{k};
	case 'extlim'
	    k = k + 1;
	    extlim = varargin{k};
	case 'suckdb'
	    k = k + 1;
	    suckdb = varargin{k};
	end
	k = k + 1;

    end

    df = fin(2) - fin(1);	% freq step

    fin0   = fin;
    snpin0 = snpin;

    kint = 1;

    while (fin0(end) * extlim < fout(end))

	% intermediate extrapolation

	fintmax    = fin0(end) * extlim;
	nfint      = round(fintmax / df);
	fint{kint} = [0:nfint]' * df;

	[snpint{kint},dnpint{kint}] = interpsnp(fin0, snpin0, fint{kint}, method, force, zerof, causal, ...
		    'maxforceg', maxge, varargin{:});

	% replace with the intermediate result
	fin0   = fint{kint};
	snpin0 = snpint{kint};
	kint   = kint + 1;
    end

    % final extrapolation

    [snpout,dnp] = interpsnp(fin0, snpin0, fout, method, force, zerof, causal, 'maxforceg', maxge, varargin{:});

    % adjust phase to the original input

    switch(force)
    case {'force','forcep'}
	[n0,n1,n2] = size(snpin);
	foutx1 = find(fout <= fin(end));
	for k1 = 1:n1
	    for k2 = 1:n2
		s1porg  = interp1(fin, snpin(:,k1,k2), fout(foutx1), 'linear', 'extrap');
		s1pout  = snpout(foutx1,k1,k2);
		s1mporg = snp2minph(fout(foutx1), s1porg) * 0;
		s1mpout = snp2minph(fout, snpout(:,k1,k2)) * 0;
		ksuck   = find(db(s1porg) < db(s1pout(1)) - suckdb, 1);
		if (length(ksuck) == 0)
		    if (fout(1) == 0)
			pderr = unwrap(angle(s1porg(2:end))) ./ (- 2*pi*fout(foutx1(2:end))) ...
			      - unwrap(angle(s1pout(2:end))) ./ (- 2*pi*fout(foutx1(2:end)));
			derr  = (1 ./ fout(foutx1(7:end))) \ (pderr(6:end) ./ fout(foutx1(7:end)));
		    else
			pderr = unwrap(angle(s1porg(1:end))) ./ (- 2*pi*fout(foutx1(1:end))) ...
			      - unwrap(angle(s1pout(1:end))) ./ (- 2*pi*fout(foutx1(1:end)));
			derr  = (1 ./ fout(foutx1(6:end))) \ (pderr(6:end) ./ fout(foutx1(6:end)));
		    end
		else
		    if (fout(1) == 0)
			pderr = unwrap(angle(s1porg(2:ksuck-1))) ./ (- 2*pi*fout(foutx1(2:ksuck-1))) ...
			      - unwrap(angle(s1pout(2:ksuck-1))) ./ (- 2*pi*fout(foutx1(2:ksuck-1)));
			derr  = (1 ./ fout(7:ksuck-1)) \ (pderr(6:ksuck-2) ./ fout(7:ksuck-1));
		    else
			pderr = unwrap(angle(s1porg(1:ksuck-1))) ./ (- 2*pi*fout(foutx1(1:ksuck-1))) ...
			      - unwrap(angle(s1pout(1:ksuck-1))) ./ (- 2*pi*fout(foutx1(1:ksuck-1)));
			derr  = (1 ./ fout(6:ksuck-1)) \ (pderr(6:ksuck-1) ./ fout(6:ksuck-1));
		    end
		end
		snpout(:,k1,k2) = snpout(:,k1,k2) .* exp(- j * fout * 2 * pi * derr);
	    end
	end
    end

end
