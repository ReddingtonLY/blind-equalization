%% read s4p
% chil_imp = snp2imp(param,2);
% figure;plot(chil_imp);
osr = 32;
method   = 'spline';
prelen   = 0.004;
dcintlen = 0.2;
suckdb   = [20 200];    % thru and XT

fs = symRate * osr;	% sampling frequency
df = 10e6;			% delta frequency
nf = round(fs/2/df);		% number of frequency points except DC
fx = [0:nf]'/nf * fs/2;		% extended frequency vector



f = [0:10e6:50e9]';
addpath snptool
chdir1 = 'E:\OneDrive - Credo Semiconductor, Inc\MATLAB project\RXequalization\test code\channel snp\bp channel\24_28_30_include_BGA_via\CaBP_BGAVia_Opt2_32dB\CaBP_BGAVia_Opt2_32dB_THRU.s4p';
[fch1, s4pch1] = loadsnp(chdir1);
s4pch1 = interpsnp(fch1, s4pch1, f);

chdir2 = 'E:\OneDrive - Credo Semiconductor, Inc\MATLAB project\RXequalization\test code\channel snp\bp channel\Cable_BKP_16dB_0p575m\Cable_BKP_16dB_0p575m_thru1.s4p';
[fch2, s4pch2] = loadsnp(chdir2);
s4pch2 = interpsnp(fch2, s4pch2, f);

chdir3 = 'E:\OneDrive - Credo Semiconductor, Inc\MATLAB project\RXequalization\test code\channel snp\bp channel\24_28_30_include_BGA_via\CaBP_BGAVia_Opt2_24dB\CaBP_BGAVia_Opt2_24dB_THRU.s4p';
[fch3, s4pch3] = loadsnp(chdir3);
s4pch3 = interpsnp(fch3, s4pch3, f);


s4p   = cats4p(s4pch1,s4pch2);%;%

%  s4p   = cats4p(s4pch1,s4pch2, s4pch3);%s4pch1;%

s4pm  = s4p2s4pm(s4p);
s21dd = s4pm(:,2,1);

% snp to impulse no xtalk
f = [0:10e6:50e9]';%fx = f;
s21ddx = extrapsnp(f, s21dd, fx, method, 'prelen', prelen, 'dcintlen', dcintlen, 'suckdb', suckdb(1));
[inp,tnp,dnp,snp,~] = snp2impl(fx, s21ddx, 'prelen', prelen, 'dcintlen', dcintlen); % convert to impulse response