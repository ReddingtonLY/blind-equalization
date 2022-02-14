function [f,snpout] = loadsnp(fn, varargin)
    snpobj = sparameters(fn);

    f     = snpobj.Frequencies;
    snpin = shiftdim(snpobj.Parameters, 2);
    	% 1st dim : frequency
	% 2nd dim : output port
	% 3rd dim : input port

    port_order  = [];
    Zref_in     = snpobj.Impedance;
    Zref_out    = snpobj.Impedance;
    Rterm       = snpobj.Impedance;
    dorenormin  = false;
    dorenormout = false;

    if (length(varargin) >= 1)
	port_order = varargin{1};	% port numbers of [IPOS INEG OPOS ONEG]
	argk = 2;
	while (argk <= length(varargin))
	    switch (varargin{argk})
	    case 'Zrefin'		% override Zref in touch stone file
		Zref_in = varargin{argk + 1};
		if (Zref_in ~= snpobj.Impedance)
		    dorenormin = true;
		end
		argk = argk + 1;
	    case 'Zrefout'
		Zref_out = varargin{argk + 1};
		if (Zref_out ~= snpobj.Impedance)
		    dorenormout = true;
		end
		argk = argk + 1;
	    case 'Rterm'
		Rterm = varargin{argk + 1};
		if (Rterm ~= snpobj.Impedance)
		    dorenormin  = true;
		    dorenormout = true;
		end
		argk = argk + 1;
	    end
	    argk = argk + 1;
	end
    end

    if (dorenormin || dorenormout)
	if (dorenormin)
	    % renormalize whole SNP from Zref_in to Rterm
	    snpterm = renormsnp(snpin, Zref_in, Rterm);
	else
	    snpterm = snpin;
	end
	% slice S4P in the convention of [IPOS INEG OPOS ONEG] = [1 3 2 4]
	s4pterm = snpterm(:, port_order([1 3 2 4]), port_order([1 3 2 4]));
	if (dorenormout)
	    % renormalize S4P from Rterm to Zref_out
	    snpout = renormsnp(s4pterm, Rterm, Zref_out);
	else
	    snpout = s4pterm;
	end
    else
	% no need to renormalize
	if (length(port_order) == 0)
	    % return whole SNP
	    snpout = snpin;
	else
	    % slice S4P in the convention of [IPOS INEG OPOS ONEG] = [1 3 2 4]
	    snpout = snpin(:, port_order([1 3 2 4]), port_order([1 3 2 4]));
	end
    end
end
