function convs32p(basename)
    [dmy,dmy] = mkdir(basename);
    disp(sprintf('loading %s.s32p\n', basename));
    s32pobj = sparameters(sprintf('%s.s32p', basename));
    disp(sprintf('saving %s/*.s4p\n', basename));
    txpos = [31 1 27 5 23  9 19 13];
    txneg = [29 3 25 7 21 11 17 15];
    rxpos = [32 2 28 6 24 10 20 14];
    rxneg = [30 4 26 8 22 12 18 16];
    f     = s32pobj.Frequencies;
    Z0    = s32pobj.Impedance;
    s32p  = s32pobj.Parameters;
    for tx = [1:8]
	for rx = [1:8]
	    ports = [txpos(tx) rxpos(rx) txneg(tx) rxneg(rx)];
	    s4p = snp2smp(s32p, Z0, ports, Z0);
	    if (tx == rx)
		fname = sprintf('%s/THRU_TX%d_RX%d.s4p', basename, tx, rx);
	    else
		fname = sprintf('%s/FEXT_TX%d_RX%d.s4p', basename, tx, rx);
	    end
	    rfwrite(s4p, f, fname);
	end
    end

end
