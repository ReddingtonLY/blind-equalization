function [s4pall, s4p] = compkg2s4p(f,p)

    % parameters
    %  p.Rd     : Tx/Rx terminator (ohm)
    %  p.Cd     : Tx/Rx device capacitor (F)
    %  p.Ls     : Tx/Rx device inductor (H)
    %  p.Cb     : Tx/Rx device bump capacitor (F)
    %  p.pkg1zp : Tx/Rx package trace 1 length (mm)
    %  p.pkg1Zc : Tx/Rx package trace 1 impedance (ohm)
    %  p.pkg2zp : Tx/Rx package trace 2 length (mm)
    %  p.pkg2Zc : Tx/Rx package trace 2 impedance (ohm)
    %  p.Cp     : Tx/Rx package capacitor (F)
    %  p.pcbC0  : Tx/Rx PCB capacitor 0 (F)
    %  p.pcbzp  : Tx/Rx PCB trace length (mm)
    %  p.pcbZc  : Tx/Rx PCB trace impedance (ohm)
    %  p.pcbC1  : Tx/Rx PCB capacitor 1 (F)

    s4p.Rd  = rdtx2snp(f, 4, p.Rd);
    s4p.Cd  = cap2snp(f, 4, p.Cd);
    s4p.Ls  = ind2snp(f, 4, p.Ls);
    s4p.Cb  = cap2snp(f, 4, p.Cb);
    s4p.TL1 = pkgck2snp(f, 4, p.pkg1zp, p.pkg1Zc);
    s4p.TL2 = pkgck2snp(f, 4, p.pkg2zp, p.pkg2Zc);
    s4p.Cp  = cap2snp(f, 4, p.Cp);

    if (p.pcbzp > 0)
	s4p.PCBC0 = cap2snp(f, 4, p.pcbC0);
	s4p.PCBTL = pcbck2snp(f, 4, p.pcbzp, p.pcbZc);
	s4p.PCBC1 = cap2snp(f, 4, p.pcbC1);
	s4p.PCB   = cats4p(s4p.PCBC0, s4p.PCBTL, s4p.PCBC1);
    else
	s4p.PCB = thru2snp(f, 4);
    end

    s4p.all = cats4p(s4p.Rd, s4p.Cd, s4p.Ls, s4p.Cb, s4p.TL1, s4p.TL2, s4p.Cp, s4p.PCB);

    s4pall = s4p.all;

end
