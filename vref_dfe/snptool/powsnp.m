function snpout = powsnp(snpin, pow)

    lnpin  = log(snpin);
    gainin = real(lnpin);		% = log(abs(snpin));
    phsin  = unwrap(imag(lnpin));	% = unwrap(angle(snpin));

    gainout = gainin * pow;
    phsout  = phsin  * pow;

    lnpout  = gainout + j * phsout;

    snpout  = exp(lnpout);

end

