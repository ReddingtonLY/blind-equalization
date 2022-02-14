%% run rx equalization without NRZ case
%clc
clear all
%close all
% setting
param.pltLvl            = 1;
param.bitRate           = 53.125*1e9*2;     %1/28.9e9;%18.82e-12 
param.osr               = 32;
param.iteTime           = 100;               %iteration times
param.mode = 'pam4';
param.eqMode = 'sep'; % com, sep
% ffe
param.ffeRange          = 6;        % scalar or vector
param.mainCursor        = 3;         % scalar or vector
% dfe
param.Ndfe              = 1;           % length of dfe
param.bmax              = 0.5;         % dfe magnitude limit

% 
param.trainingMode = 0;              % if 0, using slicer on DFE, otherwise dk_hat = training seq
param.dk_hatMode   = 1;              % if 0, use dk to update the tap weights, otherwise us dk_hat*
param.en_txfir     = 0;
param.txfir = [-0.0104    0.0532   -0.1824    0.4645   -0.2894];%[-0.0899   -0.2229    0.7351   -0.1888   -0.1140   -0.0940];
param.sampre       = 16;
% LMS, C-CMA
% c-cma_t is now for testing, pls use c-cma_t ***
[best_ffecoeff,rocoeff,dfecoeff,best_decode,y_ui_s] = rxeq(param,'lms');
best_ffecoeff = best_ffecoeff';

%% get linear fit pulse response
% settings
% param.Dp=p(1);              % equal to numPre Eq85-5: Shifted bits by D_p, denote the delay of pulse response
% param.Np=p(2);              % num rows of X1
% 
% param.Dw = p(3);
% param.Nw = p(4);            % length of ffe(3ck, 3 pre, 1 post)
% 
% param.Nv = p(5);            % 200->3ck 162.9.3.1.2 ||  13->120D 3.1.4
% param.Nb = p(6);            % 120D-8 table  Decision feedback equalizer (DFE) length
% param.eq_type = p(7);          % equalization type 0: zf, 1: mmse
% param.numPost = p(8)
% 
% OP.tx_fir_tuning = p(9);
%   Dp Np  Dw Nw Nv  Nb eq_type numPost tx_fir_tuning
p = [15 200 3  5  200 12 1       10       1];
[y_avg,P,sigma_e,X,x,c0] = get_pr(y_ui_s,best_decode,p);

figure;plot([0:1:length(P(:))-1]/32,P(:));xlabel('UI');ylabel('voltage');title('pulse response');hold on;
p_t = P(:);
stem(c0/32,p_t(c0));

%% debug
param.en_txfir = 0;param.iteTime           = 0; param.pltLvl            = 0;
% without TXFIR
p = [15 200 3  5  200 12 1       10       1];
[~,~,~,best_decode,y_ui_s] = rxeq(param,'lms');
[y_avg,P_1,sigma_e,X,x,c0] = get_pr(y_ui_s,best_decode,p);

% with TXFIR
param.en_txfir = 1;
[~,~,~,best_decode,y_ui_s] = rxeq(param,'lms');
[y_avg,P_2,sigma_e,X,x,c0] = get_pr(y_ui_s,best_decode,p);

% apply rxfir
% 162.9.3.1.1 Linear fit to the measured waveform Rm(j,i+4) MNp by 5
% matrix
param.osr = 32; % oversampling rate
if mod(param.osr,2) == 0
    m = -param.osr/2:1:param.osr/2-1;
else
    m = -(param.osr-1)/2:1:(param.osr-1)/2-1;
end
param.Np = 200; param.Dp = 15;
param.numPre  = 2;%param.mainCursor-1;
param.numPost = 8;%param.ffeRange - param.mainCursor;

P_I = zeros(param.Np*param.osr,1);P_I(15*param.osr+1:16*param.osr) = ones(param.osr,1);
[norm_tx_coeff,Rmt] = findtxcoeff(m,param,P_2(:),P_2(:));
P_3 = best_ffecoeff*Rmt';
figure;
plot([0:1:length(P_3(:))-1]/32,P_3(:));grid on


figure;
plot([0:1:length(P_1(:))-1]/32,P_1(:));hold on;
plot([0:1:length(P_2(:))-1]/32,P_2(:));
plot([0:1:length(P_3(:))-1]/32,P_3(:));
% stem(15,P_1(32,15));
% stem(16,P_1(32,16));
% stem(16.5,P_3(32*16.5));stem(17.5,P_3(32*17.5));stem(18.5,P_3(32*18.5));
xlabel('UI');ylabel('voltage');title('pulse response at RX');
legend('without TXFIR','with TXFIR','with TXFIR&RXFFE');grid on

