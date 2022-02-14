function [y_avg,P,sigma_e,X,x,cursor_i] = get_pr(y_ui,x,p)
param.symRate = 53.125e9;
param.Dp=p(1);              % Eq85-5: Shifted bits by D_p, denote the delay of pulse response
param.Np=p(2);              % num rows of X1

param.Dw = p(3);
param.Nw = p(4);            % length of ffe(3ck, 3 pre, 1 post)

param.Nv = p(5);            % 200->3ck 162.9.3.1.2 ||  13->120D 3.1.4
param.Nb = p(6);            % 120D-8 table  Decision feedback equalizer (DFE) length
param.eq_type = p(7);          % equalization type 0: zf, 1: mmse

OP.tx_fir_tuning = p(9);
OP.pl_lvl = 0;              % 1: plot equalized pulse response, and eye

param.alph = [-3 -1 1 3];   % pam4 signal level

param.numPre = p(1);           % Dp = numPre
param.numPost = p(8);


preset = 1;
[param.samplePerUI, TotalnumSym] = size(y_ui);

x = 3*x; % denormalized
if mod(param.samplePerUI,2) == 0
    m = -param.samplePerUI/2:1:param.samplePerUI/2-1;
else
    m = -(param.samplePerUI-1)/2:1:(param.samplePerUI-1)/2-1;
end

% dump garbage
%y_ui = y_ui(:,end-2001+1:end);
%x = x(end-2001+1:end);

alph = [-3 -1 1 3];
param.pattern_length = length(x);
% force repeat time to 1;
param.repeat_time = 1;
vdiff(preset) = max(y_ui(:))-min(y_ui(:));
%% cancel the DC offset
y_wo_dc = y_ui - mean(y_ui(:));



%% avg of measured waveform
% align
[RLM,P,X,x_pattern,~,m,param,waveformshift] = linearFitv2(y_wo_dc(:,1:length(x)),alph,param,x);
y_wo_dc = reshape(circshift(y_wo_dc(:),-waveformshift),param.samplePerUI,[]);

y_avg = y_wo_dc;


%% linear fit pulse response -> steady-state voltage vf and pmax
P_t = P(:);

vf(preset) = abs(sum(P_t(1:param.samplePerUI*param.Nv))/param.samplePerUI);
[vp(preset),tp] = max(P_t(:));% PMAX FOR LINEAR FIT PULSE RESPONSE
e = (P*X-y_avg);
sigma_e(preset) = sqrt(var((e(:))));cursor_i = 0;
%[cursor_i,prec,post] = findcursor(P_t,param);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [cursor_i,prec,post] = findcursor(P_t,param)
    Nosr = param.samplePerUI;
    numPre = param.numPre;
    numPost = param.numPost;

    [~, sbr_peak_i]=max(abs(P_t));
    zxi = find(diff(sign(P_t-.01*max(P_t)))>=1);
    zxi = zxi(zxi<sbr_peak_i);
    zxi = zxi(sbr_peak_i - zxi < 4*Nosr);
    if length(zxi)>1
        zxi = zxi(end);
    end
    mm_range = zxi+(0:2*Nosr);
    mm_metric = abs(P_t(mm_range-Nosr)-max(P_t(mm_range+Nosr),0));
    [~,mm_cursor_offset] = min(mm_metric);
    cursor_i = zxi+mm_cursor_offset-1;
    prec = cursor_i-[Nosr:Nosr:numPre*Nosr];
    post = cursor_i+[Nosr:Nosr:numPost*Nosr];
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RLM,P,X,x_pattern,y_ui,m,param,waveformshift]=linearFitv2(y_ui,alph,param,x_pattern)

% 162.9.3.1.1 Linear fit to the measured waveform Rm(j,i+4) MNp by 5
% matrix
if mod(param.samplePerUI,2) == 0
    m = -param.samplePerUI/2:1:param.samplePerUI/2-1;
else
    m = -(param.samplePerUI-1)/2:1:(param.samplePerUI-1)/2-1;
end

[~,PatLen] = size(y_ui);

x=x_pattern;
%% shift waveform
[waveformshift,~] = findStartIndex(y_ui,x,param);

y_ui_1 = reshape(circshift(y_ui(:),-waveformshift),param.samplePerUI,[]);


%% resample the voltage level to get normalized pattern
v = zeros(1,4);
for ii = 1:1:4
    v(ii) = mean(mean(y_ui_1(param.samplePerUI/2:param.samplePerUI/2+1,strfind(x,alph(ii))))); % sampled at middle
end
vmid = (v(1)+v(4))/2;
ES1 = (v(2)-vmid)/(v(1)-vmid);
ES2 = (v(3)-vmid)/(v(4)-vmid);
RLM = min([3*ES1,3*ES2,2-3*ES1,2-3*ES2]);
ES = (ES1+ES2)/2;

% normalized x(n) as in 3ck 162.9.3.1.1
x_norm = x;
x_norm((find(x==alph(1)))) = -1;
x_norm((find(x==alph(2)))) = -ES;
x_norm((find(x==alph(3)))) = ES;
x_norm((find(x==alph(4)))) = 1;
%*************************************************************************%
%*************************************************************************%

%% linear fit use normalized x(n) as in 3ck 162.9.3.1.1

xr = circshift(x_norm,-param.Dp);  %left shift
X = zeros(param.Np,length(x));
X(1,:) = xr;
%% X1
for ii = 1:1:param.Np-1
    X(ii+1,:) = circshift(xr,ii);
end
X(end+1,:)=ones(1,PatLen);


P = y_ui_1/X;



%% alignment start index
% Yuan@2021/1/19
% changed@2021/2/5
% y_ui: measured waveform
% x: test pattern
% os_ui: samples per UI
% N: number of samples

%function [exactIndexyuan,exactIndexyasuo, index] = findStartIndex(y_ui,x,os_ui,N)
function [ws,index] = findStartIndex(y_ui,x,param)
os_ui = param.samplePerUI;
N = param.pattern_length;
% ZERO HOLD
ideal_pulse = zeros(N*os_ui,1);
vpp_measured = [max(y_ui(:))-min(y_ui(:))];
vpp_ideal = [max(x)-min(x)];
scale = vpp_measured/vpp_ideal;
for ii=1:1:N
    ideal_pulse((ii-1)*os_ui+1:ii*os_ui)=x(ii);
end
ideal_pulse = ideal_pulse .* scale;
% cross-correlation
[coeff,lags] = xcorr(y_ui(:),ideal_pulse);
[~,id] = max(coeff);
ws = id;% waveformshift right
index = round(abs(lags(id))/os_ui)+1;
%debug
% ideal_pulseshift = circshift(ideal_pulse,(id));
% plot(ideal_pulseshift(32*1105:32*1125));hold on
% grid on
% y_t = y_ui(:);
% plot(y_t(32*1105:32*1125),'+');

%[~,exactIndexyasuo] = aligntpat(y_ui,x);


function [wave,shiftcnt] = aligntpat(wave, tpat)
wave = wave(:);
N    = length(wave);
Npn  = length(tpat);
Nosr = N / Npn;

Nsub = 4;

wavex2 = [wave; wave];
acf1 = zeros(Npn*Nsub,1);
for k = 1:Npn*Nsub
    Pmid = wavex2((k-1)*round(Nosr/Nsub)+[1:Nosr:N]);
    Pmid = Pmid - mean(Pmid);
    %	acf1(k) = abs(mean(Pmid .* tpat));
    c = corrcoef([Pmid tpat.']);
    acf1(k) = abs(c(2,1));
end

[dmy,bestk1] = max(acf1);

wvshift1 = circshift(wave, - (bestk1-1) * round(Nosr/Nsub));

tpatX = filter(ones(Nosr,1),1,upsample(tpat,Nosr));

acf2 = zeros(Nosr*2-1,1);
for k = [1:Nosr*2-1]
    Pmid = circshift(wvshift1, k - Nosr);
    Pmid = Pmid - mean(Pmid);
    %	acf2(k) = norm(Pmid .* tpatX);
    c = corrcoef([Pmid tpatX.']);
    acf2(k) = abs(c(2,1));
end

[dmy,bestk2] = max(acf2);

wvshift2 = circshift(wvshift1, bestk2 - Nosr);

wave = wvshift2;

shiftcnt = - (bestk1-1) * round(Nosr/Nsub) + bestk2 - Nosr;



%% read measurement file
% 1/18/2021 Yuan

function [y_ui,read_data]  = readMeasurement(waveform_dir)

read_data.date = 1;
fidin=fopen(waveform_dir);
row_i = 1;col_i = 1;idx = 1;
while ~feof(fidin)
    tline=fgetl(fidin);
    if isempty(tline)
        continue
    end
    
    if isempty(find(isletter(tline(1:5))))
        y_ui(row_i,col_i) = str2double(tline);
        row_i = row_i+1;
        if mod(row_i,str2double(result(15))) == 1
            row_i = 1;
            col_i = col_i+1;
        end
        continue
    else
        rst = textscan(tline,'%s','delimiter',',');
        data = strjoin(cellstr(rst{1,1}(2:end,1)));
        result{idx}=data;
        idx = idx + 1;
        
        continue
    end
end
read_data.FileFormat =              result(1);
read_data.FormatVersion =           result(2);
read_data.Instrument =              result(3);
read_data.SwVersion =               result(4);
read_data.serialNumber =            result(5);
read_data.Date =                    result(6);
read_data.sourceName =              result(7);
read_data.Points =                  str2double(result(8));
read_data.SignalTypr =              result(9);
read_data.ChannelBw =               str2double(result(10));
read_data.ChannelNoise =            str2double(result(11));
read_data.XOrg =                    str2double(result(12));
read_data.XInc =                    str2double(result(13));
read_data.SymbolRate =              str2double(result(14));
read_data.SamplesPerUI =            str2double(result(15));
read_data.PatternLength =           str2double(result(16));
read_data.X_unit =                  result(17);
read_data.Y_unit =                  result(18);

fclose(fidin);




% read prbs file
% 1/18/2021 Yuan


function [y_ui,read_data]  = readprbs(waveform_dir)

read_data.date = 1;
fidin=fopen(waveform_dir);
idxy = 1;idx = 2;
while ~feof(fidin)
    tline=fgetl(fidin);
    if isempty(tline)
        continue
    end
    
    if isempty(find(isletter(tline(1:5))))
        y_ui(idxy) = str2double(tline);
        idxy = idxy +1;
        
        continue
    else
        rst = textscan(tline,'%s','delimiter',',');
        name = regexprep(cellstr(rst{1,1}(1)),'[( ) \s /]','');
        name = name{1,1};
        if contains(tline,'File Format')
            read_data(1).file_format = name;
        else
            data = strjoin(cellstr(rst{1,1}(2:end,1)));
            read_data(idx)=struct(name, {data});
            idx = idx + 1;
        end
        continue
    end
end


fclose(fidin);

function [waveformDir, symRate,chBw,chNoise]=readWaveform(path,fname)
% default test preset is from 1 to 5b, automatic generate dir 1-5b
% the data reading is from 1 to 5b
% default - len_offset:
len_offset = 0;
defPreset = {'1'};%{'1';'2';'3';'4';'5a';'5b'};

ind = strfind(path,'/');    %linux
if isempty(ind)
    ind = strfind(path,'\');% windows
end
new_folder_wavform = [pwd path(ind(end-1):end-2) 'waveform'];
contains = dir([new_folder_wavform '/*.mat']);
%% skip read data, in order
count = 1;
for k = 1:length(contains)
    if ~strcmp(contains(k).name,[path(ind(end-1)+1:end-2) defPreset{k+count-1} '.mat'])
        puntid(count) = k;count = count+1;
    end
end
if isempty(contains) % no file read
    readlist = defPreset;
    dir_str = string(zeros(1,length(readlist)));
    for ii = 1:length(readlist)
        dir_str(ii) = strcat(path(1:end-2),readlist(ii));
    end
    symRate = zeros(1,length(dir_str));
    mkdir(new_folder_wavform);
    waveformDir = string(zeros(1,length(dir_str)));
    
    for ii = 1:1:length(dir_str)
        cdir = dir_str{ii};
        waveform_dir = sprintf('%s/%s',cdir,fname);
        file_name = sprintf('%s/%s.mat',new_folder_wavform,cdir(ind(end-1)+1:end));
        [y_ui,read_data]=readMeasurement(waveform_dir);
        
        % save to .mat
        save(file_name,'y_ui')
        symRate(ii) = read_data.SymbolRate;
        waveformDir(ii) = file_name;
    end
    
    
    
elseif ~isempty(contains)&&~(count==1) % have some file left to read
    dir_str = string(zeros(1,length(puntid)));
    for ii = 1:length(puntid)
        dir_str(ii) = strcat(path(1:end-2),defPreset(puntid(ii)));
    end
    mkdir(new_folder_wavform);
    waveformDir = string(zeros(1,length(dir_str)));
    for ii = 1:1:length(dir_str)
        cdir = dir_str{ii};
        waveform_dir = sprintf('%s/%s',cdir,fname);
        file_name = sprintf('%s/%s.mat',new_folder_wavform,cdir(ind(end-1)+1:end));
        [y_ui,~]=readMeasurement(waveform_dir);
        
        % save to .mat
        save(file_name,'y_ui');
        waveformDir(ii) = file_name;
    end
end

% only read parameter
% all data has been read but other parameter need refresh anyway
dir_str = string(zeros(1,length(defPreset)));
waveformDir = string(zeros(1,length(defPreset)));
symRate = zeros(1,length(dir_str));
chBw    = zeros(1,length(dir_str));
chNoise = zeros(1,length(dir_str));
for ii = 1:length(defPreset)
    dir_str(ii) = strcat(path(1:end-2),defPreset(ii));
end

for ii = 1:1:length(dir_str)
    cdir = dir_str{ii};
    waveform_dir = sprintf('%s/%s',cdir,fname);
    file_name = sprintf('%s/%s.mat',new_folder_wavform,cdir(ind(end-1)+1:end));
    p = readPara(waveform_dir);
    symRate(ii) = p.symbolRate;
    chBw(ii)    = p.channelBw;
    chNoise(ii) = p.channelNoise;
    waveformDir(ii) = file_name;
end



function r = readPara(path)

fidin=fopen(path);
idx = 1;
while ~feof(fidin)
    tline=fgetl(fidin);
    if isempty(tline)
        continue
    end
    
    if ~isempty(find(isletter(tline(1:5))))
        rst = textscan(tline,'%s','delimiter',',');
        data = strjoin(cellstr(rst{1,1}(2:end,1)));
        result{idx}=data;
        idx = idx + 1;
        if contains(strjoin(cellstr(rst{1,1}(1,1))),'symbol rate','IgnoreCase',true)
            r.symbolRate = str2double(data);% last one to read
            break
        elseif contains(strjoin(cellstr(rst{1,1}(1,1))),'channel bandwidth','IgnoreCase',true)
            r.channelBw = str2double(data);
            continue
        elseif contains(strjoin(cellstr(rst{1,1}(1,1))),'channel noise','IgnoreCase',true)
            r.channelNoise = str2double(data);
            continue
        end
    end
end




function [outputpath]=readTestPattern(path,name)

full_path = [path name];
%%
ind = strfind(name,'.');
new_folder = sprintf('%s/%s',pwd,name(1:ind-1));

contains = dir(new_folder);
read_flag = ~isempty(contains);
for k = 1:length(contains)
    if (contains(k).isdir)
        continue;
    end
    read_flag = strcmp(contains(k).name,[name(1:ind-1) '.mat']);
end
if ~read_flag
    [x_1,~] = readprbs(full_path);
    mkdir(new_folder);
    outputpath = [new_folder ['/' name(1:ind-1) '.mat']];
    save(outputpath,'x_1'); % default
end
outputpath = [new_folder ['/' name(1:ind-1) '.mat']];





%     %% debug
%     %% X1
%     for ii = 1:1:Np-1
%     Xd(ii+1,:) = circshift(x,ii);
%     end
%     Pa = y_avg/Xd;
%     plot(P(:),'-+');hold on
%     plot(Pa(:),'--');
%     %% debug