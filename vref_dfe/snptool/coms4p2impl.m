function [i21dd,t21dd,b21dd,d21dd,s21dd,f,ch_s21dd,rxf_s1p,ch_s21ddin,s4pmin,s4p] = coms4p2impl(fin,ch_s4p,fb,OSR,txp,rxp,varargin)

    % parameters
    %  txp.Rd, rxp.Rd         : Tx/Rx terminator (ohm)
    %  txp.Cd, rxp.Cd         : Tx/Rx device capacitor (F)
    %  txp.Ls, rxp.Ls         : Tx/Rx device inductor (H)
    %  txp.Cb, rxp.Cb         : Tx/Rx device bump capacitor (F)
    %  txp.pkg1zp, rxp,pkg1zp : Tx/Rx package trace 1 length (mm)
    %  txp.pkg1Zc, rxp,pkg1Zc : Tx/Rx package trace 1 impedance (ohm)
    %  txp.pkg2zp, rxp,pkg2zp : Tx/Rx package trace 2 length (mm)
    %  txp.pkg2Zc, rxp,pkg2Zc : Tx/Rx package trace 2 impedance (ohm)
    %  txp.Cp, rxp.Cp         : Tx/Rx package capacitor (F)
    %  txp.pcbzp, rxp,pcbzp   : Tx/Rx PCB trace length (mm)
    %  txp.pcbZc, rxp,pcbZc   : Tx/Rx PCB trace impedance (ohm)
    %  txp.c(1:5)             : Tx FIR coefficients
    %  rxp.fr                 : Rx 3dB bandwidth (BW4 filter cut-off) (Hz)
    %  rxp.gdc(1:2)           : Rx CTLE DC gains (dB)
    %  rxp.fz(1)              : Rx CTLE zero for gdc(1) = 0dB (Hz)
    %  rxp.fz(2)              : Rx CTLE zero for gdc(2) = 0dB (Hz)
    %  rxp.fp(1:3)            : Rx CTLE pole (Hz)

    s4p.txrd  = rdtx2snp(fin, 4, txp.Rd);
    s4p.txcd  = cap2snp(fin, 4, txp.Cd);
    s4p.txls  = ind2snp(fin, 4, txp.Ls);
    s4p.txcb  = cap2snp(fin, 4, txp.Cb);
    s4p.txtl1 = pkgck2snp(fin, 4, txp.pkg1zp, txp.pkg1Zc);
    s4p.txtl2 = pkgck2snp(fin, 4, txp.pkg2zp, txp.pkg2Zc);
    s4p.txcp  = cap2snp(fin, 4, txp.Cp);

    if (txp.pcbzp > 0)
	s4p.txpcb = pcbtl2snp(fin, 4, txp.pcbzp, txp.pcbZc);
    else
	s4p.txpcb = thru2snp(fin, 4);
    end

    if (rxp.pcbzp > 0)
	s4p.rxpcb = pcbtl2snp(fin, 4, rxp.pcbzp, rxp.pcbZc);
    else
	s4p.rxpcb = thru2snp(fin, 4);
    end

    s4p.rxcp  = cap2snp(fin, 4, rxp.Cp);
    s4p.rxtl2 = pkgck2snp(fin, 4, rxp.pkg2zp, rxp.pkg2Zc);
    s4p.rxtl1 = pkgck2snp(fin, 4, rxp.pkg1zp, rxp.pkg1Zc);
    s4p.rxcb  = cap2snp(fin, 4, rxp.Cb);
    s4p.rxls  = ind2snp(fin, 4, rxp.Ls);
    s4p.rxcd  = cap2snp(fin, 4, rxp.Cd);
    s4p.rxrd  = rdtx2snp(fin, 4, rxp.Rd);

    s4pin = cats4p(s4p.txrd, s4p.txcd, s4p.txls, s4p.txcb, s4p.txtl1, s4p.txtl2, s4p.txcp, s4p.txpcb, ch_s4p, ...
		   s4p.rxpcb, s4p.rxcp, s4p.rxtl2, s4p.rxtl1, s4p.rxcb, s4p.rxls, s4p.rxcd, s4p.rxrd);

    fs = fb * OSR;		% sampling frequency = symbol rate * oversampling ratio
    df = round(diff(fin(1:2)));	% frequency pitch
    nf = round(fs/2/df);	% number of frequency grid - 1
    f  = [0:nf]'/nf*fs/2;	% frequency grid for interpolation (Hz=1/sec)

    s4pmin = s4p2s4pm(s4pin);					% convert S4P to mixed-mode S4P

    ch_s21ddin = s4pmin(:,2,1);					% S21dd with package

    ch_s21dd = interpsnp(fin, ch_s21ddin, f, varargin{:});	% extrapolate channel S21dd

    if (rxp.fr > 0)

	% apply RX filter

	bw4_s1p = buttern2s1p(f, rxp.fr, 4);			% Rx noise filter (BW4 filter)

	w  = 2 * pi * f;		% angular frequency (rad/s)
	s  = w * j;			% imaginary angular frequency (rad/s)

	ctle_s1p =  (10^(rxp.gdc(1)/20) + s / (2 * pi * rxp.fz(1))) ...
		 .* (10^(rxp.gdc(2)/20) + s / (2 * pi * rxp.fz(2))) ...
		 ./ (1 + s / (2 * pi * rxp.fp(1)))                  ...
		 ./ (1 + s / (2 * pi * rxp.fp(2)))                  ...
		 ./ (1 + s / (2 * pi * rxp.fp(3)));			% Rx CTLE

	rxf_s1p = bw4_s1p .* ctle_s1p;				% concatenate Rx noise filter and CTLE

	s21dd = ch_s21dd .* rxf_s1p;				% concatenate channel and Rx filter

    else

	s21dd = ch_s21dd;					% channel only

    end

    [i21dd,t21dd,d21dd] = snp2impl(f, s21dd);		% convert S21dd to impulse response

    ui = 1/fb;
    
    u21dd  = cumsum(i21dd);	% step response
    b21dd0 = u21dd - interp1(t21dd, u21dd, t21dd - ui, 'linear', 0);	% SBR

    b21dd  = zeros(length(b21dd0), 1);

    for k = [1:length(txp.c)]
	b21dd_tap = interp1(t21dd, b21dd0, t21dd - ui * (k-1), 'linear', 0);	% SBR of each FIR tap
	b21dd     = b21dd + b21dd_tap * txp.c(k);
    end

end
