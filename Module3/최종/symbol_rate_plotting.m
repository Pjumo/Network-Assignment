clearvars; clc;
%--------Simulation parameters----------------
nSym=10^4; %Number of OFDM Symbols to transmit
EbN0dB = 0:2:20; % bit to noise ratio
MOD_TYPE='MPSK'; %modulation type - 'MPSK' or 'MQAM'
% M --> PSK: [2 4 8 16 32 64]
% M --> QAM: [4 16 64 256]
M=4; %choose modulation order for the chosen MOD_TYPE
N=64; %FFT size or total number of subcarriers (used + unused) 64
Ncp= 16; %number of symbols in the cyclic prefix
L=10; %Number of taps for the frequency selective channel model
%--------Derived Parameters--------------------
k=log2(M); %number of bits per modulated symbol
EsN0dB = 10*log10(k*N/(N+Ncp))+EbN0dB; %account for overheads
errors1= zeros(1,length(EsN0dB)); %to store symbol errors RAYLEIGH
errors2= zeros(1,length(EsN0dB)); %to store symbol errors RICIAN
K_dB = 3;   %K factor in DB
K = 10.^(K_dB/10); %K factor in linear scale
g1 = sqrt(K/(2*(K+1))); g2 = sqrt(1/(2*(K+1)));

for i=1:length(EsN0dB)%Monte Carlo Simulation
    for j=1:nSym
    %-----------------Transmitter--------------------
    d=ceil(M.*rand(1,N));%uniform distributed random syms from 1:M
    [X,ref]=modulation_mapper(MOD_TYPE,M,d);
    x= ifft(X,N);% IDFT
    s = add_cyclic_prefix(x,Ncp); %add CP
    %-------------- Channel ----------------
    h_Ray =1/sqrt(2)*(randn(1,L)+1i*randn(1,L)); %RAYLEIGH Fading 
    h_Ric = (g2*randn(1,L)+g1)+1i*(g2*randn(1,L)+g1);% RICIAN Fading
    H_Ray = fft(h_Ray,N);
    H_Ric = fft(h_Ric,N);
    hs_Ray = conv(h_Ray,s);%filter the OFDM sym through freq. sel. channel
    hs_Ric = conv(h_Ric,s);
    r_Ray = add_awgn_noise(hs_Ray,EsN0dB(i));%add AWGN noise r = s+n
    r_Ric = add_awgn_noise(hs_Ric,EsN0dB(i));%add AWGN noise r = s+n
    %-----------------Receiver----------------------
    y_Ray = remove_cyclic_prefix(r_Ray,Ncp,N);%remove CP
    y_Ric = remove_cyclic_prefix(r_Ric,Ncp,N);%remove CP
    Y_Ray = fft(y_Ray,N);%DFT
    Y_Ric = fft(y_Ric,N);%DFT
    
    V_Ray = Y_Ray./H_Ray;%equalization
    V_Ric = Y_Ric./H_Ric;
    [~,dcap_ray]= iqOptDetector(V_Ray,ref);%demapper using IQ detector
    [~,dcap_ric]= iqOptDetector(V_Ric,ref);%demapper using IQ detector
    %----------------Error counter------------------
    numErrors_Ray=sum(d~=dcap_ray);%Count number of symbol errors
    numErrors_Ric=sum(d~=dcap_ric);%Count number of symbol errors
    errors1(i)=errors1(i)+numErrors_Ray;%accumulate symbol errors
    errors2(i)=errors2(i)+numErrors_Ric;%accumulate symbol errors
    end
end

% PLOT SER for Rayleigh and Rician Fading Channels

SER_sim_RAY = errors1/(nSym*N);%symbol error from simulation RAYLEIGH fading
SER_theory_RAY=ser_rayleigh(EbN0dB,MOD_TYPE,M); %theoretical SER for Rayleigh flat-fading
figure();
semilogy(EbN0dB,SER_sim_RAY,'ko');hold on;
semilogy(EbN0dB,SER_theory_RAY,'k-');grid on;
title(['Performance of ',num2str(M),'-', MOD_TYPE,' OFDM over RAYLEIGH channel']);
xlabel('Eb/N0 (dB)');ylabel('Symbol Error Rate');
legend('simulated','theoretical');


SER_theory_RIC = ser_rician(EbN0dB,K_dB,MOD_TYPE,M);%theoretical SER Rician
SER_sim_RIC = errors2/(nSym*N);%symbol error from simulation RIC fading
figure();
semilogy(EbN0dB,SER_sim_RIC,'ko');hold on;
semilogy(EbN0dB,SER_theory_RIC,'k-');grid on;
title(['Performance of ',num2str(M),'-', MOD_TYPE,' OFDM over RICIAN channel']);
xlabel('Eb/N0 (dB)');ylabel('Symbol Error Rate');
legend('simulated','theoretical');

% ========================= SUPPORTING FUNCTIONS (no need to edit) =========================%


function [ser] = ser_rician(EbN0dB,K_dB,MOD_TYPE,M)
%Compute Theoretical Symbol Error rates for MPSK or MQAM modulations
%EbN0dB - list of SNR per bit points
%K_dB - K factor for Rician fading in dB
%MOD_TYPE - 'MPSK' or 'MQAM'
%M - Modulation level for the chosen modulation
% - For MPSK M can be any power of 2
% - For MQAM M must be even power of 2 (square QAM only)
gamma_b = 10.^(EbN0dB/10); %SNR per bit in linear scale
gamma_s = log2(M)*gamma_b; %SNR per symbol in linear scale
K=10^(K_dB/10); %K factor in linear scale
switch lower(MOD_TYPE)
    case {'mpsk','psk'}
        ser = zeros(size(gamma_s));
        for i=1:length(gamma_s) %for each SNR point
            g = sin(pi/M).^2;
            fun = @(x) ((1+K)*sin(x).^2)/((1+K)*sin(x).^2+g*gamma_s(i)).*exp(-K*g*gamma_s(i)./((1+K)*sin(x).^2+g*gamma_s(i))); %MGF
            ser(i) = (1/pi)*integral(fun,0,pi*(M-1)/M);
        end
    case {'mqam','qam'}
        ser = zeros(size(gamma_s));
        for i=1:length(gamma_s) %for each SNR point
            g = 1.5/(M-1);
            fun = @(x) ((1+K)*sin(x).^2)/((1+K)*sin(x).^2+g*gamma_s(i)).*exp(-K*g*gamma_s(i)./((1+K)*sin(x).^2+g*gamma_s(i))); %MGF
            ser(i) = 4/pi*(1-1/sqrt(M))*integral(fun,0,pi/2)-4/pi*(1-1/sqrt(M))^2*integral(fun,0,pi/4);
        end
end
end

function [dCap]=demodulate(MOD_TYPE,M,r,COHERENCE)
%Wrapper functin to call various digital demodulation techniques
% MOD_TYPE - 'PSK','QAM','PAM','FSK'
% M - modulation order, For BPSK M=2, QPSK M=4, 256-QAM M=256 etc..,
% COHERENCE - only applicable if FSK modulation is chosen
% - 'coherent' for coherent MFSK
% - 'noncoherent' for coherent MFSK
% r - received modulated symbols
% dCap - demodulated information symbols
%Construct the reference constellation for the selected mod. type
switch lower(MOD_TYPE)
    case 'bpsk'
        dCap= mpsk_detector(2,r); %M=2
    case 'mpsk'
        dCap= mpsk_detector(M,r);
    case 'mqam'
        dCap= mqam_detector(M,r);
    case 'pam'
        dCap= mpam_detector(M,r);
    case 'fsk'
        dCap= mfsk_detector(M,r,COHERENCE);
    otherwise
        error('Invalid Modulation specified');
end
end

function [dCap]= mpsk_detector(M,r)
%Function to detect MPSK modulated symbols
%[dCap]= mpsk_detector(M,r) detects the received MPSK signal points
%points - 'r'. M is the modulation level of MPSK
ref_i= 1/sqrt(2)*cos(((1:1:M)-1)/M*2*pi);
ref_q= 1/sqrt(2)*sin(((1:1:M)-1)/M*2*pi);
ref = ref_i+1i*ref_q; %reference constellation for MPSK
[~,dCap]= iqOptDetector(r,ref); %IQ detection
end

function [X,ref]=modulation_mapper(MOD_TYPE,M,d)
%Modulation mapper for OFDM transmitter
% MOD_TYPE - 'MPSK' or 'MQAM' modulation
% M - modulation order, For BPSK M=2, QPSK M=4, 256-QAM M=256 etc..,
% d - data symbols to be modulated drawn from the set {1,2,...,M}
%returns
% X - modulated symbols
% ref -ideal constellation points that could be used by IQ detector
if strcmpi(MOD_TYPE,'MPSK')
    [X,ref]=mpsk_modulator(M,d);%MPSK modulation
else
    if strcmpi(MOD_TYPE,'MQAM')
        [X,ref]=mqam_modulator(M,d);%MQAM modulation
    else
        error('Invalid Modulation specified');
    end
end
end

function [s,ref]=mqam_modulator(M,d)
%Function to MQAM modulate the vector of data symbols - d
%[s,ref]=mqam_modulator(M,d) modulates the symbols defined by the vector d
% using MQAM modulation, where M specifies order of M-QAM modulation and
% vector d contains symbols whose values range 1:M. The output s is modulated
% output and ref represents reference constellation that can be used in demod
if(((M~=1) && ~mod(floor(log2(M)),2))==0) %M not a even power of 2
    error('Only Square MQAM supported. M must be even power of 2');
end
ref=constructQAM(M); %construct reference constellation
s=ref(d); %map information symbols to modulated symbols
end

function [ref,varargout]= constructQAM(M)
%Function to construct gray codes symbol constellation for M-QAM
% [ref]=constructQAM(M) - returns the ideal signaling points (ref) in a
% symmetric rectangular M-QAM constellation, where M is the level of QAM
% modulation. The returned constellation points are arranged such that the
% index of the points are arranged in a Gray-coded manner. When plotted,
% indices of constellation points will differ by 1 bit.
%
% [ref,I,Q]=constructQAM(M) - returns the ideal signaling points (ref) along
% with the IQ components breakup in a symmetric rectangular M-QAM constellation,
% where M is the level of QAM modulation. The returned constellation points are
% arranged such that the index of the points are arranged in a Gray-coded manner.
n=0:1:M-1; %Sequential address from 0 to M-1 (1xM dimension)
%------Addresses in Kmap - Gray code walk---------------
a=dec2gray(n); %Convert linear addresses to gray code
N=sqrt(M); %Dimension of K-Map - N x N matrix
a=reshape(a,N,N).'; %NxN gray coded matrix
evenRows=2:2:size(a,1); %identify alternate rows
a(evenRows,:)=fliplr(a(evenRows,:));%Flip rows - KMap representation
nGray=reshape(a.',1,M); %reshape to 1xM - Gray code walk on KMap
%Construction of ideal M-QAM constellation from sqrt(M)-PAM
D=sqrt(M); %Dimension of PAM constellation
x=floor(nGray/D);
y=mod(nGray,D);
Ax=2*(x+1)-1-D; %PAM Amplitudes 2m-1-D - real axis
Ay=2*(y+1)-1-D; %PAM Amplitudes 2m-1-D - imag axis
ref=Ax+1i*Ay; %assign complex numbers to reference
if nargout==2 %if only one variable argument is given
    varargout{1}=Ax; %Real part (I)
elseif nargout==3 %if two variable arguments are given
    varargout{1}=Ax; %Real part (I)
    varargout{2}=Ay; %Imaginary part (Q)
end
end

function [grayCoded]=dec2gray(decInput)
%convert decimal to Gray code representation
%example: x= [0 1 2 3 4 5 6 7] %decimal
%y=dec2gray(x)
%returns y = [ 0 1 3 2 6 7 5 4] %Gray coded
[rows,cols]=size(decInput);
grayCoded=zeros(rows,cols);
for i=1:rows
    for j=1:cols
        grayCoded(i,j)=bitxor(bitshift(decInput(i,j),-1),decInput(i,j));
    end
end
end

function [s,ref]=mpsk_modulator(M,d)
%Function to MPSK modulate the vector of data symbols - d
%[s,ref]=mpsk_modulator(M,d) modulates the symbols defined by the
%vector d using MPSK modulation, where M specifies the order of
%M-PSK modulation and the vector d contains symbols whose values
%in the range 1:M. The output s is the modulated output and ref
%represents the reference constellation that can be used in demod
ref_i= 1/sqrt(2)*cos(((1:1:M)-1)/M*2*pi);
ref_q= 1/sqrt(2)*sin(((1:1:M)-1)/M*2*pi);
ref = ref_i+1i*ref_q;
s = ref(d); %M-PSK Mapping
end

function s = add_cyclic_prefix(x,Ncp)
%function to add cyclic prefix to the generated OFDM symbol x that
%is generated at the output of the IDFT block
% x - ofdm symbol without CP (output of IDFT block)
% Ncp-num. of samples at x's end that will copied to its beginning
% s - returns the cyclic prefixed OFDM symbol
s = [x(end-Ncp+1:end) x]; %Cyclic prefixed OFDM symbol
end

function [r,n,N0] = add_awgn_noise(s,SNRdB,L)
%Function to add AWGN to the given signal
%[r,n,N0]= add_awgn_noise(s,SNRdB) adds AWGN noise vector to signal
%'s' to generate a %resulting signal vector 'r' of specified SNR
%in dB. It also returns the noise vector 'n' that is added to the
%signal 's' and the spectral density N0 of noise added
%
%[r,n,N0]= add_awgn_noise(s,SNRdB,L) adds AWGN noise vector to
%signal 's' to generate a resulting signal vector 'r' of specified
%SNR in dB. The parameter 'L' specifies the oversampling ratio used
%in the system (for waveform simulation). It also returns the noise
%vector 'n' that is added to the signal 's' and the spectral
%density N0 of noise added
s_temp=s;
if iscolumn(s), s=s.'; end %to return the result in same dim as 's'
gamma = 10^(SNRdB/10); %SNR to linear scale
if nargin==2, L=1; end %if third argument is not given, set it to 1
if isvector(s)
    P=L*sum(abs(s).^2)/length(s);%Actual power in the vector
else %for multi-dimensional signals like MFSK
    P=L*sum(sum(abs(s).^2))/length(s); %if s is a matrix [MxN]
end
N0=P/gamma; %Find the noise spectral density
if(isreal(s))
    n = sqrt(N0/2)*randn(size(s));%computed noise
else
    n = sqrt(N0/2)*(randn(size(s))+1i*randn(size(s)));%computed noise
end
r = s + n; %received signal
if iscolumn(s_temp), r=r.'; end%return r in original format as s
end

function y = remove_cyclic_prefix(r,Ncp,N)
%function to remove cyclic prefix from the received OFDM symbol r
% r - received ofdm symbol with CP
% Ncp - num. of samples at beginning of r that need to be removed
% N - number of samples in a single OFDM symbol
% y - returns the OFDM symbol without cyclic prefix
y=r(Ncp+1:N+Ncp);%cut from index Ncp+1 to N+Ncp
end

function [dCap]= mqam_detector(M,r)
%Function to detect MQAM modulated symbols
%[dCap]= mqam_detector(M,r) detects the received MQAM signal points
%points - 'r'. M is the modulation level of MQAM
if(((M~=1) && ~mod(floor(log2(M)),2))==0) %M not a even power of 2
    error('Only Square MQAM supported. M must be even power of 2');
end
ref=constructQAM(M); %reference constellation for MQAM
[~,dCap]= iqOptDetector(r,ref); %IQ detection
end

function [idealPoints,indices]= iqOptDetector(received,ref)
%Optimum Detector for 2-dim. signals (MQAM,MPSK,MPAM) in IQ Plane
%received - vector of form I+jQ
%ref - reference constellation of form I+jQ
%Note: BPSK are one dim. modulations. The same function can be
%applied for these modulations since quadrature is zero (Q=0).
x=[real(received); imag(received)]';%received vec. in cartesian form
y=[real(ref); imag(ref)]';%reference vec. in cartesian form
[idealPoints,indices]= minEuclideanDistance(x,y);
end

function [idealPoints,indices]= minEuclideanDistance(x,y)
%function to compute the pairwise minimum Distance between two
%vectors x and y in p-dimensional signal space and select the
%vectors in y that provides the minimum distances.
% x - a matrix of size mxp
% y - a matrix of size nxp. This acts as a reference against
% which each point in x is compared.
% idealPoints - contain the decoded vector
% indices - indices of the ideal points in reference matrix y
[m,p1] = size(x);[n,p2] = size(y);
if p1~=p2
    error('Dimension Mismatch: x and y must have same dimension')
end
X = sum(x.*x,2);
Y = sum(y.*y,2)';
d = X(:,ones(1,n)) + Y(ones(1,m),:) - 2*x*y';%Squared Euclidean Dist.
[~,indices]=min(d,[],2); %Find the minimum value along DIM=2
idealPoints=y(indices,:);
indices=indices.';
end

function [ser] = ser_rayleigh(EbN0dB,MOD_TYPE,M)
%Compute Theoretical Symbol Error rates for MPSK or MQAM modulations
%EbN0dB - list of SNR per bit points
%MOD_TYPE - 'MPSK' or 'MQAM'
%M - Modulation level for the chosen modulation
% - For MPSK M can be any power of 2
% - For MQAM M must be even power of 2 (square QAM only)
gamma_b = 10.^(EbN0dB/10); %SNR per bit in linear scale
gamma_s = log2(M)*gamma_b; %SNR per symbol in linear scale
switch lower(MOD_TYPE)
    case {'bpsk'}
        ser = 0.5*(1-sqrt(gamma_b/(1+gamma_b)));
    case {'mpsk','psk'}
        ser = zeros(size(gamma_s));
        for i=1:length(gamma_s) %for each SNR point
            g = sin(pi/M).^2;
            fun = @(x) 1./(1+(g.*gamma_s(i)./(sin(x).^2))); %MGF
            ser(i) = (1/pi)*integral(fun,0,pi*(M-1)/M);
        end
    case {'mqam','qam'}
        ser = zeros(size(gamma_s));
        for i=1:length(gamma_s) %for each SNR point
            g = 1.5/(M-1);
            fun = @(x) 1./(1+(g.*gamma_s(i)./(sin(x).^2)));%MGF
            ser(i) = 4/pi*(1-1/sqrt(M))*integral(fun,0,pi/2)-4/pi*(1-1/sqrt(M))^2*integral(fun,0,pi/4);
        end
    case {'mpam','pam'}
        ser = zeros(size(gamma_s));
        for i=1:length(gamma_s) %for each SNR point
            g = 3/(M^2-1);
            fun = @(x) 1./(1+(g.*gamma_s(i)./(sin(x).^2)));%MGF
            ser(i) = 2*(M-1)/(M*pi)*integral(fun,0,pi/2);
        end
end
end