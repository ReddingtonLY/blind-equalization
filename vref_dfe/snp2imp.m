
%% Load S4P and transform to impulse response
function [chil_imp] = snp2imp(param,chidx)

addpath snptool

fprintf('Constructing channel\n');

% parameters to convert S-parameter to impulse response

%method   = 'linear';
method   = 'spline';
%method   = 'makima';
prelen   = 0.004;
dcintlen = 0.2;
suckdb   = [20 200];    % thru and XT

fs = param.bitRate * 2;	% sampling frequency
df = 10e6;			% delta frequency
nf = round(fs/2/df);		% number of frequency points except DC
fx = [0:nf]'/nf * fs/2;		% extended frequency vector

switch (chidx)

case 1
    % 1: PCIe6_64G_HPE_pad_to_pad_channels_2020_03_05/S4P/Channel1_TL85Ohm_32dB

    f = [0:10e6:60e9]';	% this is same as fch of channel below

    chdir = 'D:/cygwin64/chdata/public/PCIE/PCIe6_64G_HPE_pad_to_pad_channels_2020_03_05/S4P/Channel1_TL85Ohm_32dB/';
    s4pfn = {'THRU.s4p', 'FEXT1.s4p', 'FEXT2.s4p', 'NEXT1.s4p', 'NEXT2.s4p', 'NEXT3.s4p'};
    Nfext = chnfext(chidx);
    Nnext = chnnext(chidx);

    PCB4dBZc  = 100;
    PCB4dBzp  = 150;
    s4pPCB4dB = pcbck2snp(f, 4, PCB4dBzp, PCB4dBZc);

    % S-parameter files of termination network

    txtmrs4pfn = './projectfiles/viceroy/S4P/io_network_tx.s4p';
    rxtmrs4pfn = './projectfiles/viceroy/S4P/io_network_rx.s4p';

    [ftxtmr, s4ptxtmr] = loadsnp(txtmrs4pfn);	% [0:5e6:100e9]
    s4ptx = interpsnp(ftxtmr, s4ptxtmr, f);

    [frxtmr, s4prxtmr0] = loadsnp(rxtmrs4pfn);	% [0:5e6:100e9]
    %frxtmr = frxtmr(2:end);
    %s4prxtmr = renormsnp(s4prxtmr(2:end,:,:), 50, [50 50e8 50 50e8]);
    s4prxtmr = renormsnp(s4prxtmr0, 50, [50 50e8 50 50e8]);
    s4prxtmr(:,2,1) = s4prxtmr(:,2,1) * 1e4;
    s4prxtmr(:,4,3) = s4prxtmr(:,4,3) * 1e4;
    s4prxtmr(:,1,2) = s4prxtmr(:,1,2) * 1e4;
    s4prxtmr(:,3,4) = s4prxtmr(:,3,4) * 1e4;
    s4prx = interpsnp(frxtmr, s4prxtmr, f);

    for k = [1:length(s4pfn)]
	[fch, s4pch] = loadsnp([chdir s4pfn{k}]);	% fch is same as f
	if (k <= 1 + Nfext)
	    % THRU or FEXT channel
	    s4p   = cats4p(s4pPCB4dB, s4ptx, s4pch, s4prx);	% concatenate 4dB PCB loss and TX and RX terminators
	    s4pm  = s4p2s4pm(s4p);
	    s21dd = s4pm(:,2,1);
	    if (k == 1)
		% THRU channel
		s21ddx = extrapsnp(f, s21dd, fx, method, 'prelen', prelen, 'dcintlen', dcintlen, 'suckdb', suckdb(1));
	    else
		% FEXT channel
		s21ddx = extrapsnp(f, s21dd, fx, method, 'prelen', prelen, 'dcintlen', dcintlen, 'suckdb', suckdb(2));
	    end
	else
	    % NEXT channel
	    s4p   = cats4p(s4ptx, s4pch, s4prx);	% concatenate TX and RX terminators
	    s4pm  = s4p2s4pm(s4p);
	    s21dd = s4pm(:,2,1);
	    s21ddx = extrapsnp(f, s21dd, fx, method, 'prelen', prelen, 'dcintlen', dcintlen, 'suckdb', suckdb(2));
	end

	impl{k} = snp2impl(fx, s21ddx, 'prelen', prelen, 'dcintlen', dcintlen); % convert to impulse response
    end

case {2,3}
    % 2: C2M 2 inch + COM 12mm package
    % 3: C2M 9 inch + COM 12mm package

    f = [0:10e6:60e9]';

    % COM device and package model

    pkgp.Rd     =  50;		%  p.Rd     : Tx/Rx terminator (ohm)
    pkgp.Cd     = 120e-15;	%  p.Cd     : Tx/Rx device capacitor (F)
    pkgp.Ls     = 120e-12;	%  p.Ls     : Tx/Rx device inductor (H)
    pkgp.Cb     =  30e-15;	%  p.Cb     : Tx/Rx device bump capacitor (F)
    pkgp.pkg1zp =  12;		%  p.pkg1zp : Tx/Rx package trace 1 length (mm)
    pkgp.pkg1Zc =  87.5;	%  p.pkg1Zc : Tx/Rx package trace 1 impedance (ohm)
    pkgp.pkg2zp =   1.8;	%  p.pkg2zp : Tx/Rx package trace 2 length (mm)
    pkgp.pkg2Zc =  92.5;	%  p.pkg2Zc : Tx/Rx package trace 2 impedance (ohm)
    pkgp.Cp     =  87e-15;	%  p.Cp     : Tx/Rx package capacitor (F)
    pkgp.pcbC0  =   0e-15;	%  p.pcbC0  : Tx/Rx PCB capacitor 0 (F)
    pkgp.pcbzp  =   0;		%  p.pcbzp  : Tx/Rx PCB trace length (mm)
    pkgp.pcbZc  = 100;		%  p.pcbZc  : Tx/Rx PCB trace impedance (ohm)
    pkgp.pcbC1  =   0e-15;	%  p.pcbC1  : Tx/Rx PCB capacitor 1 (F)

    s4ptx = compkg2s4p(f, pkgp);
    s4prx = s4ptx(:,[2 1 4 3],[2 1 4 3]);

    setpathC2M;

    switch (chidx)
    case 2
	c2mchidx = 114; % c2mch114 : JLchan5/lim_3ck_adhoc_02_073119/Channel5a_Smaller_Pad_2inch_trace/
    case 3
	c2mchidx = 117; % c2mch117 : JLchan5/lim_3ck_adhoc_02_073119/Channel5d_Smaller_Pad_9inch_trace/
    end

    chdir = c2mchpath{c2mchidx}.dir;
    s4pfn{1} = c2mchpath{c2mchidx}.thru;
    s4pfn{2} = c2mchpath{c2mchidx}.fext{1};
    s4pfn{3} = c2mchpath{c2mchidx}.fext{2};
    s4pfn{4} = c2mchpath{c2mchidx}.fext{3};
    s4pfn{5} = c2mchpath{c2mchidx}.fext{4};
    s4pfn{6} = c2mchpath{c2mchidx}.fext{5};
    s4pfn{7} = c2mchpath{c2mchidx}.next{1};
    s4pfn{8} = c2mchpath{c2mchidx}.next{2};
    s4pfn{9} = c2mchpath{c2mchidx}.next{3};
    Nfext = 5;%chnfext(chidx);
    Nnext = 3;%chnnext(chidx);

    for k = [1:length(s4pfn)]
	[fch, orgs4pch] = loadsnp([chdir s4pfn{k}]);	% fch is [50e6:10e6:50e9] (2") or [10e6:10e6:50e9] (9")
	s4pch = interpsnp(fch, orgs4pch, f);		% align frequency
	s4p   = cats4p(s4ptx, s4pch, s4prx);		% concatenate TX and RX package models
	s4pm  = s4p2s4pm(s4p);
	s21dd = s4pm(:,2,1);
	if (k == 1)
	    % THRU channel
	    s21ddx = extrapsnp(f, s21dd, fx, method, 'prelen', prelen, 'dcintlen', dcintlen, 'suckdb', suckdb(1));
	else
	    % FEXT or NEXT channel
	    s21ddx = extrapsnp(f, s21dd, fx, method, 'prelen', prelen, 'dcintlen', dcintlen, 'suckdb', suckdb(2));
	end
	impl{k} = snp2impl(fx, s21ddx, 'prelen', prelen, 'dcintlen', dcintlen); % convert to impulse response
    end
end

%% Channel
% plot(fx*1.e-9,20*log10(abs(s21ddx)));xlabel('Ghz');ylabel('dB');grid on


chil_imp = impl{1};
chil_imp = chil_imp(1:round(length(chil_imp)/2));

%%% FEXT and NEXT

% clear orgchfext;
% for Fi=1:1:Nfext
%     orgchfext{Fi} = impl{1+Fi};
% end
% clear orgchnext;
% for Ni=1:1:Nnext
%     orgchnext{Ni} = impl{1+Nfext+Ni};
% end

