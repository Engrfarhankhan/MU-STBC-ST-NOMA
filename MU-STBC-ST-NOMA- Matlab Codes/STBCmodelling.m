%% MU-STBC-ST

warning off;
clear all;
clc

Nbps=1; 
M=2^Nbps;  % Modulation order=2/4/6 for QPSK/16QAM/64QAM
EbN0=0:5:25; 
Eb_N0=-15:5:15;    % EbN0
N_iter=30000;   % no of packets   % Number of iterations for each EbN0
Nframe=1;         % no. of sybmols in on % Number of symbols per frame
sigPow=0;         % Signal power initialization
norms=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];     % BPSK 4-QAM 16-QAM

%% Bob
PowerdB=[0 -8 -17 -21 -25]; % Channel tap power profile 'dB'
Delay=[0 3 5 6 8];          % Channel delay 'sample'
Power=10.^(PowerdB/10);     % Channel tap power profile 'linear scale'
Power=Power/(sum(Power));
Ntap=length(PowerdB);       % Chanel tap number
Lch=Delay(end)+1;           %Channel length


Nfft=64;           % FFT size
Ng=16; % Guard interval length
Nsym=Nfft+Ng;      % Symbol duration
Nvc=0;        % Nvc=0: no virtual carrier
Nused=Nfft ; %-Nvc;

NgType=1; % NgType=1/2 for cyclic prefix/zero padding
if NgType==1, nt='CP';  elseif NgType==2, nt='ZP';   end
% EbN0Lin = 10.^(EbN0/10);
% channel_H = [];channel_H2 = [];
% p = 0.02872; q_new = 1.5;
% p = 0.02872; q = 1.4;
% 
% errors= zeros(1,length(EbN0));
% errors_new= zeros(1,length(EbN0));
%%


for k=1:length(EbN0)
      BERu1=0;
       BERu2=0;
     BERu1s=0;
       BERu2s=0;
     
     
     TERrx1=0;
     TERrx2=0;
        TERrx1s=0;
     TERrx2s=0;
%      TotalTER=0;
     
     bittot=0;
     bitTot=0;

     
     for m=1:N_iter
      
      %Antenna 1,
      %Tx1 to Rx1       
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h11=zeros(1,Lch); 
       h11(Delay+1)=channel; % cir: channel impulse response
        H11=fft([h11 zeros(1,Nfft-Lch)]); % Channel frequency response
        H11 = diag(H11);
       
      %Tx1 to Rx2       
        channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h12=zeros(1,Lch); 
       h12(Delay+1)=channel; % cir: channel impulse response
%        H12A2=fft(h12,64);
        H12=diag(fft([h12 zeros(1,Nfft-Lch)])); % Channel frequency response
        
     
 
 
 
      %Antenna 2, 
      %Tx2 to Rx1
       channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h21=zeros(1,Lch); 
       h21(Delay+1)=channel; % cir: channel impulse response
       H21=diag(fft([h21 zeros(1,Nfft-Lch)])); % Channel frequency response
       
      %Tx2 to Rx2
        channel=(randn(1,Ntap)+j*randn(1,Ntap)).*sqrt(Power/2);%%%%%
       h22=zeros(1,Lch); 
       h22(Delay+1)=channel; % cir: channel impulse response

       H22=diag(fft([h22 zeros(1,Nfft-Lch)])); % Channel frequency response 64x64 diagonal matrix
      
        
     
%% 
% * Generating Signal for each antenna, modulating, and adding CP*      

       % X1 Data
      X= randi([0,M-1], Nused*Nframe,1)'; 
      Xmod= transpose(qammod(X,M,'gray')/norms(Nbps)); % 64x1  1,-1	 
      U1=Xmod;	 
      % X2 Data
      X2= randi([0,M-1], Nused*Nframe,1)'; 
      Xmod2= transpose(qammod(X2,M,'gray')/norms(Nbps));% 64x1 1,-1
      U2=Xmod2;
%       UT=XmodT;
 
%% 
% *Calculating Auxiliary Signals a1, a2*

     
%       a2 = ((H11A1 + H12A2).*((H21A1 + H22A2).*P1+(H23A3 + H24A4).*P2).*U1 - (H21A1 + H22A2).*((H11A1 + H12A2).*P1+(H13A3 + H14A4).*P2).*U2)./((H13A3 + H14A4).*(H21A1 + H22A2)-(H11A1 + H12A2).*(H23A3 + H24A4));
%       a1 = (-((H21A1 + H22A2).*P1+(H23A3 + H24A4).*P2).*U1-(H23A3 + H24A4).*a2)./(H21A1 + H22A2);
% 
      w = H12 - H11.*inv(H21).*H22;
      a2 = inv(w).*((H11 + H11.*inv(H21).*H22).*U1 - (H11 + H12).*U2);
      a1 = -(eye(64) + inv(H21).*H22).*U1-(inv(H21).*H22.*a2);
      
      
      w = H12 - H11.*inv(H21).*H22;
      a4 = inv(w).*((-H11 + H11.*inv(H21).*H22).*conj(U1) - (-H11 + H12).*conj(U2));
      a3 = -(-eye(64) + inv(H21).*H22).*conj(U1)-(inv(H21).*H22.*a4);
%% 
% *Superpositioning*

      %Time Slot 1
      %Transmit Antenna 1,

%       Tx1=(U1 + U2) + diag(a1);
      Tx1=(U1 + U2)+ diag(a1);
      %Transmit Antenna 2,
      
      Tx2=(U1 + U2)+ diag(a2);
      
      %Time Slot 2
       %Transmit Antenna 1,

      Tx1s=(-conj(U2)-conj(U1)) + diag(a3);
   
      %Transmit Antenna 2,
      
      Tx2s=(conj(U1)+conj(U2))+ diag(a4);
     
%% 
% *For PAPR*

       x1np= ifft(U1);
       x_GI1np= guard_interval(Ng,Nfft,NgType,x1np);       
       x2np= ifft(U2);
       x_GI2np= guard_interval(Ng,Nfft,NgType,x2np);
%% 

      % this from antenna 1   
      x1= ifft(Tx1)/norm(ifft(Tx1),2); %% This is how i am normalizing it :step1
      x_GI1= guard_interval(Ng,Nfft,NgType,x1); % add guard bands  
      % this from antenna 2
      x2= ifft(Tx2)/norm(ifft(Tx2),2);  
      x_GI2= guard_interval(Ng,Nfft,NgType,x2);
      
      % this from antenna 1   
      x1s= ifft(Tx1s)/norm(ifft(Tx1s),2); %% This is how i am normalizing it :step1
      x_GI1s= guard_interval(Ng,Nfft,NgType,x1s); % add guard bands  
      % this from antenna 2
      x2s= ifft(Tx2s)/norm(ifft(Tx2s),2);  
      x_GI2s= guard_interval(Ng,Nfft,NgType,x2s);
      
      
       % convolution for the signal from antenna 1 in first time slot, u1
      y11= conv(x_GI1,h11); 
      y21= conv(x_GI1,h21); 
   % convolution for the signal from antenna 1 in second time slot, u1
      y11s= conv(x_GI1s,h11); 
      y21s= conv(x_GI1s,h21); 
%% 

      % convolution for the signal from antenna 2 in first time slot, u2
      y12= conv(x_GI2,h12); 
      y22= conv(x_GI2,h22); 
      
      % convolution for the signal from antenna 2 in second time slot, u2
      y12s= conv(x_GI2s,h12); 
      y22s= conv(x_GI2s,h22); 
        
         %%%%%%% FOR PAPR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         if k==1
         PAPRtx1(1,m)=(max((abs(x_GI1)).^2)/mean((abs(x_GI1)).^2));
         PAPRtx2(1,m)=(max((abs(x_GI2)).^2)/mean((abs(x_GI2)).^2));
         PAPRtx1s(1,m)=(max((abs(x_GI1s)).^2)/mean((abs(x_GI1s)).^2));
         PAPRtx2s(1,m)=(max((abs(x_GI2s)).^2)/mean((abs(x_GI2s)).^2));
         PAPRtx3(1,m)=max((abs(x_GI1np)).^2)/mean((abs(x_GI1np)).^2);
         PAPRtx4(1,m)=max((abs(x_GI2np)).^2)/mean((abs(x_GI2np)).^2);
         %pause();
         end
         
         %%%%%%% FOR PAPR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
               
         snr = EbN0(k)+10*log10(Nbps*(Nused/Nfft)); % SNR vs. Eb/N0
       %noise during first time slot
         EEs=sum(x1).^2/length(x1);
         EEs2 = sum(x2).^2/length(x2);
         
         EEs_average = (EEs + EEs2)*0.5;
         sigma_noise  =sqrt(0.5*EEs_average*((30.^(-snr/5))));
         
       %noise during second time slot
         EE1s=sum(x1s).^2/length(x1s);
         EE2s = sum(x2s).^2/length(x2s);
         
         EE_averages = (EE1s + EE2s)*0.5;
         sigma_noises  =sqrt(0.5*EE_averages*((30.^(-snr/5))));
      %For noise in first time slot  
      y1n=y11+y12;
      y2n=y21+y22;
      %For noise in second time slot  
      y1ns=y11s+y12s;
      y2ns=y21s+y22s;
%       Signals Received During First Time slot (T1)
      y1 = norm(ifft(Tx1),2)*(y11) + norm(ifft(Tx2),2)*(y12) + sigma_noise.*(randn(size(y1n))+j*randn(size(y1n)));
   
      y2 = norm(ifft(Tx1),2)*(y21) + norm(ifft(Tx2),2)*(y22) + sigma_noise.*(randn(size(y2n))+j*randn(size(y2n)));

%      Signals Received During Second Time slot (T2)
      y1s = norm(ifft(Tx1s),2)*(y11s) + norm(ifft(Tx2s),2)*(y12s) + sigma_noises.*(randn(size(y1ns))+j*randn(size(y1ns)));
   
      y2s = norm(ifft(Tx1s),2)*(y21s) + norm(ifft(Tx2s),2)*(y22s) + sigma_noises.*(randn(size(y2ns))+j*randn(size(y2ns)));
         
      %During T1
      rGI1 = remove_GI(Ng,Nsym,NgType,y1); %this is time domain 
      rGI2 = remove_GI(Ng,Nsym,NgType,y2); %this is time domain 
      %During T2
      rGI1s = remove_GI(Ng,Nsym,NgType,y1s); %this is time domain 
      rGI2s = remove_GI(Ng,Nsym,NgType,y2s); %this is time domain 

      % Applying fft During T1 and Equalizing
      YU1= fft(rGI1).*diag(inv(H11 + H12));
      YU2= fft(rGI2).*diag(inv(H21 + H22));   

%       % Applying fft During T2 and Equalizing
      YU1s= fft(rGI1s).*diag(inv(-H11 + H12));
      YU2s= fft(rGI2s).*diag(inv(-H21 + H22)); 

%% 
% *Combining*

     Y1=(YU1+YU1s);
    Y2=(YU2+YU2s);

      %% Demodulation
      UU1i=qamdemod(Y1*norms(Nbps),M,'gray');
      UU2i=qamdemod(Y2*norms(Nbps),M,'gray');
      
         
%% 
% *Calculating BER and PER*

         % BER during T1
          BERu1= BERu1+sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps)));
             BERU1p(k,m)=sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps))); %for PER
          BERu2= BERu2+sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps)));
             BERU2p(k,m)=sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps))); %for PER


%% 
% *Calculating Throughput*

        % TER during T1
        TERrx1=TERrx1+(Nfft*Nbps)-sum(sum(de2bi(UU1i,Nbps)~=de2bi(X,Nbps))); %Rx1
        TERrx2=TERrx2+(Nfft*Nbps)-sum(sum(de2bi(UU2i,Nbps)~=de2bi(X2,Nbps))); %Rx2

        bittot=bittot+Nfft*Nbps;
        bitTot=bitTot+2*Nfft*Nbps;
     end


      % during T1
      BERu11(k)=BERu1/bittot;
      BERu21(k)=BERu2/bittot;
      AvagTBER(k)=(BERu1+BERu2)/bitTot;
%        % during T2
%       BERu11s(k)=BERu1s/bittot;
%       BERu21s(k)=BERu2s/bittot;
%       AvagTBERs(k)=(BERu1s+BERu2s)/bitTot;
      
%     Thorughput
      resour=N_iter*Nfft;
 
      TERrx11(k)=TERrx1/(resour);
      TERrx21(k)=TERrx2/resour;
      TotalTER(k)=(TERrx1+TERrx2)/resour;

end

N = 2*10^7; % number of bits or symbols
Eb_N0_dB = [0:25]; % multiple Eb/N0 values
nRx = 2;
for ii = 1:length(Eb_N0_dB)

    % Transmitter
    ip = rand(1,N)>0.5; % generating 0,1 with equal probability
    s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0

    % Alamouti STBC 
    sCode = 1/sqrt(2)*kron(reshape(s,2,N/2),ones(1,2)) ;

    % channel
    h = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)]; % Rayleigh channel
    n = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)]; % white gaussian noise, 0dB variance

    y = zeros(nRx,N);
    yMod = zeros(nRx*2,N);
    hMod = zeros(nRx*2,N);
    for kk = 1:nRx

        hMod = kron(reshape(h(kk,:),2,N/2),ones(1,2)); % repeating the same channel for two symbols    
        hMod = kron(reshape(h(kk,:),2,N/2),ones(1,2));
        temp = hMod;
        hMod(1,[2:2:end]) = conj(temp(2,[2:2:end])); 
        hMod(2,[2:2:end]) = -conj(temp(1,[2:2:end]));

        % Channel and noise Noise addition
        y(kk,:) = sum(hMod.*sCode,1) + 10^(-Eb_N0_dB(ii)/20)*n(kk,:);

        % Receiver
        yMod([2*kk-1:2*kk],:) = kron(reshape(y(kk,:),2,N/2),ones(1,2));
    
        % forming the equalization matrix
        hEq([2*kk-1:2*kk],:) = hMod;
        hEq(2*kk-1,[1:2:end]) = conj(hEq(2*kk-1,[1:2:end]));
        hEq(2*kk,  [2:2:end]) = conj(hEq(2*kk,  [2:2:end]));

    end

    % equalization 
    hEqPower = sum(hEq.*conj(hEq),1);
    yHat = sum(hEq.*yMod,1)./hEqPower; % [h1*y1 + h2y2*, h2*y1 -h1y2*, ... ]
    yHat(2:2:end) = conj(yHat(2:2:end));

    % receiver - hard decision decoding
    ipHat = real(yHat)>0;

    % counting the errors
    nErr(ii) = size(find([ip- ipHat]),2);
    tErr(ii) = N-size(find([ip- ipHat]),2);

end

simBer = nErr/N; % simulated ber
simTer = tErr/N; % simulated ter


EbN0Lin = 10.^(EbN0/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p));
pAlamouti = 1/2 - 1/2*(1+2./EbN0Lin).^(-1/2);
theoryBerAlamouti_nTx2_nRx1 = pAlamouti.^2.*(1+2*(1-pAlamouti)); 
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));
%% Plots

%BER
figure(1);
semilogy(EbN0,BERu11,'*r-','LineWidth',2);
hold on;
semilogy(EbN0,BERu21,'ob-','LineWidth',2);
hold on;
semilogy(EbN0,theoryBer,'*g-','LineWidth',2);   hold on;
semilogy(EbN0,theoryBerAlamouti_nTx2_nRx1,'c+-','LineWidth',2);hold on;
semilogy(Eb_N0_dB,simBer,'kd-','LineWidth',2); hold off;
title('BER Vs EbN0 for MU-STBC-ST with OFDM');
xlabel('SNR (dB)');
ylabel('Bit Error Rate(BER)');
grid on;
legend('UA-BER','UB-BER','Conv. OFDM-BER', 'Theory STBC (nTx=2, nRx=1)','Conv. STBC(nTx=2, nRx=2)');


%%
% %%TER
figure(7);
semilogy(EbN0,TERrx11,'or--','LineWidth',2);hold on; %Rx1
semilogy(EbN0,TERrx21,'*b--','LineWidth',2);hold on; %Rx2
% semilogy(Eb_N0,TERrx11s,'og--','LineWidth',2);hold on; %Rx1
% semilogy(Eb_N0,TERrx21s,'*c--','LineWidth',2);hold on; %Rx2
semilogy(EbN0,TotalTER,'mx-','LineWidth',2); hold on; %Rx1 & %Rx2
semilogy(Eb_N0_dB,simTer,'kd-','LineWidth',2); hold off;
title('TER Vs EbN0 for MU-STBC-ST with OFDM');
xlabel('SNR (dB)');
ylabel('Throughput (TER)');
grid on;
legend('UA-TER','UB-TER','TER-GAIN (UA & UB)','STBC-TER');
%%
% %PAPR
PAPRtx1dB=10*log10(PAPRtx1); % in dB
[F1,X1] = ecdf(PAPRtx1dB); %CDF

PAPRtx2dB=10*log10(PAPRtx2); % in dB
[F2,X2] = ecdf(PAPRtx2dB); %CDF

PAPRtx1dBs=10*log10(PAPRtx1s); % in dB
[F1s,X1s] = ecdf(PAPRtx1dBs); %CDF

PAPRtx2dBs=10*log10(PAPRtx2s); % in dB
[F2s,X2s] = ecdf(PAPRtx2dBs); %CDF

PAPRtx3dB=10*log10(PAPRtx3); % in dB
[F3,X3] = ecdf(PAPRtx3dB); %CDF

PAPRtx4dB=10*log10(PAPRtx4); % in dB
[F4,X4] = ecdf(PAPRtx4dB); %CDF
figure(8)
semilogy(X1,1-F1,'or--','LineWidth',2); hold on; % 1-F1 just to draw the CCDF
semilogy(X2,1-F2,'*b--','LineWidth',2); hold on;
semilogy(X1s,1-F1s,'*g-','LineWidth',2); hold on; % 1-F1 just to draw the CCDF
semilogy(X2s,1-F2s,'*m-','LineWidth',2); hold on;
semilogy(X3,1-F3,'kp-','LineWidth',2); hold on; % 1-F1 just to draw the CCDF
semilogy(X4,1-F4,'cd-','LineWidth',2); 
title('PAPR Vs CCDF for MU-STBC-ST with OFDM');
xlabel('PAPR');
ylabel('CCDF');
grid on;
legend('MU-STBC-ST-1','MU-STBC-ST-2','MU-STBC-ST-3','MU-STBC-ST-4','Conv. OFDM-1','Conv. OFDM-2');