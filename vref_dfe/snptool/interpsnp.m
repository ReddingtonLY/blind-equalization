function [snpout,dnp,s1pgerr] = interpsnp(fin,snpin,fout,varargin)

    % interpolate S parameter in different frequency

    method = 'makima';	% method of interpolation to step response in time domain
    forceg = true;	% enforce given S parameter gain in freq domain
    forcep = true;	% enforce given S parameter phase in freq domain
    zerof  = false;	% enforce zero for >x2 of given frequency
    causal = false;	% enforce causality in time-domain
    maxge  = 10;	% max gain enforcement in dB

    k = 1;
    while (k <= length(varargin))
	opt = varargin{k};
	switch(opt)
	case {'spline','linear','makima'}
	    method = opt;
	case {'force','forceg'}
	    forceg = true;
	    forcep = false;
	case 'forcep'
	    forceg = true;
	    forcep = true;
	case 'noforce'
	    forceg = false;
	    forcep = false;
	case 'zerof'
	    zerof = true;
	case 'nozerof'
	    zerof = false;
	case 'causal'
	    causal = true;
	case 'nocausal'
	    causal = false;
	case 'maxforceg'
	    k = k + 1;
	    maxge = varargin{k};
	end
	k = k + 1;

    end

    % convert to impulse response

    [inpin, tnpin, dnp] = snp2impl(fin, snpin, varargin{:});

    % convert to step response

    unpin = cumsum(inpin);

    % sampling time

    tsi = 1 / 2 / fin(end);			% time step of snpin
    tso = 1 / 2 / fout(end);			% time step of snpout
    no  = round(fout(end)/diff(fout(1:2)));
    to  = [0:2*no-1]' * tso;

    [n0,n1,n2] = size(snpin);

    tnpout = zeros(2*no,n1,n2);
    unpout = zeros(2*no,n1,n2);
    inpout = zeros(2*no,n1,n2);
    if (fout(1) == 0)
	snpout = zeros(no+1,n1,n2);
    else
	snpout = zeros(no,n1,n2);
    end

    foutx1 = find(fout <= fin(end));
    foutx2 = foutx1(end) * 2 - foutx1;
    extra  = foutx2(1) - length(fout);
    if (extra > 0)
	foutx1 = foutx1(extra+1:end);
	foutx2 = foutx2(extra+1:end);
    end
    foutx3 = find(fout > fin(end) * 2);

    for k1 = 1:n1
	for k2 = 1:n2
	    t1pin = [dnp(k1,k2) - tsi; tnpin(:,k1,k2)];
	    u1pin = [unpin(1,k1,k2);   unpin(:,k1,k2)];
	    t1pout = to + dnp(k1,k2);
	    u1pout = interp1(t1pin, u1pin, t1pout, method, u1pin(end));
	    if (causal)
		if (abs(snpin(1,k1,k2)) > 0.5)
		    % apply delay if DC gain > -6dB
		    if (fin(1) == 0)
			pd1pin = unwrap(angle(snpin(2:end,k1,k2))) ./ (- 2 * pi * fin(2:end));
		    else
			pd1pin = unwrap(angle(snpin(1:end,k1,k2))) ./ (- 2 * pi * fin(1:end));
		    end
		    fd1pin = min(pd1pin);
		else
		    fd1pin = 0;
		end
		knegout = find(t1pout < fd1pin, 1, 'last');
		if (- (t1pout(knegout) - fd1pin) >= t1pout(knegout + 1) - fd1pin)
		    u1pout(knegout) = u1pout(knegout - 1);
		end
	    end
	    i1pout = diff([0; u1pout]);
	    s1psym = fft(i1pout);
	    if (fout(1) == 0)
		s1pout = s1psym(1:no+1) .* exp(- j * fout * 2 * pi * dnp(k1,k2));
	    else
		s1pout = s1psym(2:no+1) .* exp(- j * fout * 2 * pi * dnp(k1,k2));
	    end
	    if (forceg)
		s1porg  = interp1(fin, snpin(:,k1,k2), fout(1:foutx1(end)), 'linear', 'extrap');
		s1pgerr = abs(s1porg) ./ abs(s1pout(1:foutx1(end)));
		s1pgerr(find(s1pout(1:foutx1(end)) == 0)) = 1;
		s1pgerr = min(s1pgerr, 10^(maxge/20));
		s1pout(1:foutx1(end))   = s1pout(1:foutx1(end))   .* s1pgerr;
		s1pout(foutx2(1:end-1)) = s1pout(foutx2(1:end-1)) .* s1pgerr(foutx1(1:end-1));
		if (forcep)
		    ksuck = find(db(s1porg(1:foutx1(end))) < db(s1porg(1)) - 20, 1);
		    if (length(ksuck) == 0)
			if (fout(1) == 0)
			    s1ppderr = unwrap(angle(s1porg(2:foutx1(end)))) ./ (- 2 * pi * fout(2:foutx1(end))) ...
				     - unwrap(angle(s1pout(2:foutx1(end)))) ./ (- 2 * pi * fout(2:foutx1(end)));
			else
			    s1ppderr = unwrap(angle(s1porg(1:foutx1(end)))) ./ (- 2 * pi * fout(1:foutx1(end))) ...
				     - unwrap(angle(s1pout(1:foutx1(end)))) ./ (- 2 * pi * fout(1:foutx1(end)));
			end
		    else
			if (fout(1) == 0)
			    s1ppderr = unwrap(angle(s1porg(2:ksuck-1))) ./ (- 2 * pi * fout(2:ksuck-1)) ...
				     - unwrap(angle(s1pout(2:ksuck-1))) ./ (- 2 * pi * fout(2:ksuck-1));
			else
			    s1ppderr = unwrap(angle(s1porg(1:ksuck-1))) ./ (- 2 * pi * fout(1:ksuck-1)) ...
				     - unwrap(angle(s1pout(1:ksuck-1))) ./ (- 2 * pi * fout(1:ksuck-1));
			end
		    end
		    s1pderr = mean(s1ppderr);
		    s1pout  = s1pout .* exp(- j * fout * 2 * pi * s1pderr);
		end
	    end
	    if (zerof)
		s1pout(foutx3) = s1pout(foutx3) * 1e-100;
	    end
	    tnpout(:,k1,k2) = t1pout;
	    unpout(:,k1,k2) = u1pout;
	    inpout(:,k1,k2) = i1pout;
	    snpout(:,k1,k2) = s1pout;
	end
    end
end
