
%%Preprocessing
%Load EEG and MEG data

EEG = zef_EEG.measurements;
MEG = zef_MEG.measurements;

%Transformation
EEG_Z = zscore(EEG); %student normalization
EEG_L = log(EEG);
MEG_Z = zscore(MEG);
MEG_L = log(MEG);

%Baseline correction of ERP analysis
ERP_BE(:,1) = mean(EEG_Z(:,1:30),2);
ERP_BM(:,1) = mean(MEG_Z(:,1:30),2);

for i = 2:331
    ERP_BE(:,i) = mean(EEG_Z(:,i:i+29),2);
    EEG_B(:,i) = EEG(:,i)-ERP_BE(:,i);
    PvalueE(:,i) = ttest(ERP_BE(:,i),EEG_B(:,i))
    
    ERP_BM(:,i) = mean(MEG_Z(:,i:i+29),2);
    MEG_B(:,i) = MEG(:,i)-ERP_BM(:,i);
    PvalueB(:,i) = ttest(ERP_BM(:,i),MEG_B(:,i))
end

%SEP
k = size(EEG_B,1);
m = mean(mean(EEG_B,2));
% the Nth channel
LogPower = log(abs((fft(EEG_B(n,:)))).^2);
GMFP_E = sqrt(sum(((EEG_B(n,:)-m).^2),1)/k); %Global mean field potential
figure;
subplot(2,1,1)
plot(1:360,EEG_B(n,:),'b--',1:360,GMFP_E,'r');
xlabel('Time');
ylabel('Amplitude')
legend('Log power of ERP corrected baseline','GMFP');
title('SEPs for EEG');
subplot(2,1,2)
plot(1:331,PvalueE)
title('Pvalue for EEG');


% figure,
% subplot(2,1,1)
% plot(1:360,EEG(2,1:360),'b',1:360,EEG_T(2,:),'r--')
% title('ERP baseline for EEG');
% subplot(2,1,2)
% plot(1:360,Pvalue)
% title('Pvalue for EEG');
% figure;
% k = size(MEG_T,1);
% m = mean(MEG_T,2);
% GMFP_M = sqrt(sum(((MEG_T-m).^2),1)/k);
% plot(1:360,MEG_T','b',1:360,GMFP_M,'r');
% xlabel('Time');
% ylabel('Amplitude')
% legend('Measurements','GMFP');
% title('SEPs for MEG');



%%
X(1,:) =  EEG_T(5,31:390);
X(2,:) =  EEG_T(33,31:390);
X(3,:) =  EEG_T(18,31:390);

F = fft(EEG_T(5,100:150));
for j =1:3
    for i = 1:2:51
    m_F(j,i) = mean(F(j,i:i+1),2);
    [~,p(j,i),~,~] = ttest(F(j,i),F(j,i+1));
    [~,m_p(j,i),~,~] = ttest(m_F(j,i),m_F(j,i+1));
    end
    figure,
    subplot(2,1,1)
    plot(1:53,log(F(j,:)));
    hold on
    plot(1:51,log(m_F(j,:)),'--');
    hold on
    subplot(2,1,2)
    plot(1:53,log(p(j,:)));
    hold on
    plot(1:51,log(m_F(j,:)),'--');
    hold off
end
%% time-frequency transformation(Hanning window = “raised cosine”)

M = 31;         % Window length
n = 2*M+1;
N = 75;         % FFT length (zero padding around a factor of 2)
fs = 1200;
fpass = 20;
fstop = 250;
%deltf = 3.1/M;
fc = (fstop-fpass)/2/fs; %fc is the midpoint of our transition band from 20Hz to 250Hz
%%Gamma
X=zef_MEG.measurements;
FFT_M1 = fft_hanning(X,M,n,N,fs,fpass,fstop,fc)
FFT_M2 = fft_slide(X,M,n,fs,fpass,fstop,fc,100,150)
%%
% Compute Hanning window:
X = fft(EEG_T);
pspectrum(X);
w = [2*fc*sinc(2*pi*(1:n)) zeros(1,N-n)] %hanning window
X_W = zeros(72,N+1) %initialte the zero pad out
FFT_M = zeros();
X_W(1,:) = [X(1,1:length(w)+1)] .* fft([2*fc,w]);% fft with convolution   
M_X(1,:) = mean(X_W(1,:),2); %average the slide windows
FFT_M(:,1) = mean(M_X(1,:));
FFT_T(:,1) = M_X(1,1);
for i = 2:size(X,1)-1
    for j = 1:360-N-1
    X_W(i,j:j+N) = X(i,1:length(w)+1) .* fft([2*fc,w]);% fft with convolution   
    M_X(i,j) = mean(X_W(i,j:j+N),2); %average the slide windows
    
    end
    FFT_M(:,i) = mean(M_X(i,:));
    FFT_T(:,i) = FFT_M(:,i)+FFT_M(:,i+1);
end

% Display real part of windowed signal and the Hanning window:
subplot(3,1,1)
plot(1:N,w(1,:),'-');
title('Hanning Window and Windowed, Zero-Padded');
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(3,1,2)
plot(1:N,real(FFT_M(1,:)));
title('Fast Fourier (Real Part)');
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(3,1,3)
plot(1:N,(abs(FFT_M(1,:))).^2,'ro');
plot on 
plot(1:N,(abs(FFT_M(2,:))).^2,'go');
plot on 
plot(1:N,(abs(FFT_M(3,:))).^2,'bo');
title('Power');
xlabel('Time (samples)'); ylabel('Spectrum');

%wvtool(hann(n));

%%Gamma band stimuli 
for i =1:9
figure,
subplot(3,3,i)
area(DTFX2{i});
xlim('[10 50]')
title('Monkey data recorded about saccade taske (stratium -> LPFT)')
figure,
subplot(3,3,i)
area(PDC(i));
xlim('[0 96]')
title('Subject data recorded about fingure stimulation (3a -> 1)')
end


