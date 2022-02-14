function [ERL, t, PTDR, Grr, Gloss] = s11dd2erl(f, s11dd, fb, L, Tfx, p)

% calculate ERL from S11dd

% f		Frequency vector
% s11dd		Diff-to-diff one-port s-parameter
% fb		Baudrate
% L		Number of signal levels
% Tfx		Twice the propagation delay with the test fixture
% p		Parameter structure

% parameters in structure p

% p.tr		Transition time
% p.fr		Receiver 3dB bandwidth
% p.N		Length of the reflection signal
% p.M		Number of samples per unit interval
% p.Nbx		Equalizer length associated with reflection signal
% p.betax	Incremental available signal loss factor
% p.rhox	Permitted reflection from transmission line external to DUT
% p.DER0	Target detector error ratio

% internal parameters

    halfNbin = 1000;
    prelen   = 0.040;
    dcintlen = 0.1;
    %method   = 'linear';
    method   = 'spline';
    %method   = 'makima';
    extlim   = 1.4;

    % extrapolate S11dd to high frequency

    fs  = fb * p.M;
    df  = diff(f(1:2));
    nf  = round(fs/2/df);
    fx  = [0:nf]'/nf*fs/2;	% extended frequency grid

    s11ddx = extrapsnp(f, s11dd, fx, method, 'extlim', extlim, 'prelen', prelen, 'dcintlen', dcintlen);

    Hr  = buttern2s1p(fx, p.fr, 4);					% (93A-20)
    Ht  = gauss2s1p(fx, p.tr);						% (93A-46)

    H11 = Ht .* s11ddx .* Hr;						% (93A-58)

    Tb  = 1 / fb;
    X   = Tb * sin(pi * fx * Tb) ./ (pi * fx * Tb);			% (93A-23) with At = 1
    X(find(fx == 0)) = Tb;	% avoid NaN at DC

    [PTDR, t, dly] = snp2impl(fx, X .* H11 / Tb);			% (93A-59)

    PTDR = PTDR * p.M;		% multiply for the oversampling ratio

    x1  = find((t >= Tfx) .* (t <  Tfx + (p.Nbx + 1) / fb));
    x2  = find(              (t >= Tfx + (p.Nbx + 1) / fb));
    t1  = t(x1);

    % (93A-61)

    Grr     = zeros(length(t), 1);
    Grr(x1) = p.rhox * (1 + p.rhox) * exp(- ((t1 - Tfx) * fb - (p.Nbx+1)).^2 / (p.Nbx+1)^2);
    Grr(x2) = 1;

    % (93A-62)

    Gloss     = zeros(length(t), 1);
    Gloss(x1) = 10 .^ (p.betax/fb * ((t1 - Tfx) * fb - (p.Nbx+1)) / 20);
    Gloss(x2) = 1;

    Reff = PTDR .* Grr .* Gloss;					% (93A-60)

    hlinear = Reff(find(t >= Tfx));
    hall = reshape(hlinear(1:p.M*p.N), p.M, p.N)';			% (93A-63)

    sigma = sum(hall.^2);						% (93A-64)

    % 93A.5.4 x-quantile of the reflection distribution

    [dmy, m] = max(sigma);
    h = hall(:,m);		% h(n) is column m of hall, where m maximizes sigma

    ymax = sum(abs(h));
    y    = [-halfNbin:halfNbin]'/halfNbin * ymax;
    pdf  = [zeros(halfNbin,1); 1; zeros(halfNbin,1)];

    for n = [1:p.N]
	pdf_n = zeros(halfNbin*2+1, 1);
	for l = [0:L-1]
	    pdf_n = pdf_n + DiracDelta(y, (2*l/(L-1)-1)*h(n)) / L;	% (93A-39)
	end
	pdf = conv(pdf, pdf_n, 'same');					% (93A-40)
    end

    cdf = cumsum(pdf);							% (93A-37)

    % 93A.5.5 ERL

    x = [1; 1 + find(diff(cdf) ~= 0)];	% avoid indices where cdf does not change

    ERL = - 20 * log10(- interp1(cdf(x), y(x), p.DER0, 'linear'));

end

% PDF of Dirac delta function on discrete axis
% (peak is split to two nearest grids)

function pdf = DiracDelta(yaxis, dy)
    pdf         = zeros(length(yaxis), 1);
    idx_hi      = find(yaxis > dy, 1);
    idx_lo      = idx_hi - 1;
    pitch       = yaxis(idx_hi) - yaxis(idx_lo);
    w_lo        = (yaxis(idx_hi) - dy) / pitch;
    w_hi        = 1 - w_lo;
    pdf(idx_lo) = w_lo;
    pdf(idx_hi) = w_hi;
end

