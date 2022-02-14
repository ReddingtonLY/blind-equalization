function savesnp(snp, f, fn, varargin)

    Zref = 50;

    argk = 1;
    while (argk <= length(varargin))
	switch (varargin{argk})
	case 'Zref'
	    Zref = varargin{argk + 1};
	    argk = argk + 1;
	end
	argk = argk + 1;
    end

    if (size(snp,2) == 1)
	snpout = shiftdim(snp, -2);
    else
	snpout = shiftdim(snp, 1);
    end

    if exist(fn, 'file')
	delete(fn);
    end

    rfwrite(sparameters(snpout, f, Zref), fn);

end
