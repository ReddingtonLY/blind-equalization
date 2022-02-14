
% 
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

chdir = '/nfshome/redliu/Documents/MATLAB/COM/BP_Channel/24_28_30_include_BGA_via/CaBP_BGAVia_Opt2_24dB/CaBP_BGAVia_Opt2_24dB_THRU.s4p';



f = [0:10e6:50e9]';

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

[fch, orgs4pch] = loadsnp(chdir);	% fch is [50e6:10e6:50e9] (2") or [10e6:10e6:50e9] (9")
s4pch = interpsnp(fch, orgs4pch, f);		% align frequency
s4p   = cats4p(s4ptx, s4pch, s4prx);		% concatenate TX and RX package models
s4pm  = s4p2s4pm(s4p);
s21dd = s4pm(:,2,1);

% extraplation
s21ddx = extrapsnp(f, s21dd, fx, method, 'prelen', prelen, 'dcintlen', dcintlen, 'suckdb', suckdb(1));


