function[TF,DF_num,DF_denom,DI,ka,fvec]= FFTset_v4(filename,speaker,soundspeed,a)
%Function which takes the FFT of a complete data set, and calculated
%Direcitivty Factor and Directivity Index.

%Input:
%filename - file of the dataset to be treated
%speaker - either value 1 or 2, depending on which speaker was utilized
%soundspeed - value of the velocity of sound in the medium, in [m/s]
%a - radius; used in calculating the ka values

%Output:
%TF - transfer function
%DF_num - numerator of Directivity Factor
%DF_denom - denominator of Directivity Factor
%DI - Directivity Index
%ka - frequency axis converted to ka numbers
%fvec - frequency axis in Hz

%Written by Isabel Berg, 14.10.21

%% Loading .mat file with data for measurement series
data = filename;

%% FFT
%Without edge reflection:
w_start = 250;
if speaker == 1 % Making end point for window optional:
    %w_end = 550;
    w_end = 600;
else
    %w_end = 1000;
    w_end = 600;
end
nfft = 16384;
fs = 48000;
fvec = fs/nfft*[0:nfft/2-1]; 

%Making half Hanning window:
n_end = 100;
end_win = hanning(n_end*2);
end_win(1:n_end) = [];
n_full = w_end - w_start + 1;
full_win = ones(n_full,1);
full_win((end - n_end + 1):end) = end_win(:);

TF_set = fft(data(w_start:w_end,:).*full_win,nfft);

TF_set = TF_set(1:nfft/2,:);
TF = TF_set;

%% Directivity factor and directivity index
%Since abs(TF_set) is proportional to the sound pressure amplitude for each
%frequency value of the measurements, a matrix of pressure amplitudes can
%be built: (angle in one dimension, frequency in the other)
P_matrix = TF_set(:,1:90).';

c = soundspeed;
lambda = c./fvec; %Wave length of incoming plane wave
wavenum = (2*pi)./lambda; %Wave number
ka = wavenum.*a;

%Looping through all frequencies:
DF_num = abs((P_matrix(1,:)).^2); %Indexed 1 for angle 0
DF_denom = zeros(1,length(DF_num));
for i = 1:length(fvec)
    y = sum(abs((P_matrix(:,i)).^2));
    DF_denom(1,i) = y;
end
DF_denom = DF_denom/size(P_matrix,1);

%Calculating Directivity Factor:
DF = DF_num./DF_denom;
DF_val = DF;

%Calculating Directivity Index:
 DI_val = 10*log10(DF_val);
 DI = DI_val;


