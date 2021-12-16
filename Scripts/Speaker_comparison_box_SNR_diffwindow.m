%Creates plot of the directivity factor, but with the numerator and the
%denominator plotted separately. Includes the noise floor x2. Utilizes the
%function FFTset_v4.m, as well as importing measurement, noise, and C data
%files. Written for the box scattering object.

%Original written 31.08.21 by Isabel Berg
%Updated 19.10.21 to use FFTset_v4 (including C data and radius), and make
%the code more efficient. Also changed ka number to kl number.

clear all
close all
%% REQUIRES THE FUNCTION "FFTset_v4.m"

%% REQUIRES LOCALLY STORED MEASUREMENT, NOISE, AND C DATA FILES (TOT: 24)

%% Opening data sets, using function for running FFT and calculating DF
dataset = cell(1,20);
%Loading all data sets (PLACEMENTS TO BE CHANGED WHEN RUN ELSEWHERE):
numfiles = 20;
for k = 1:numfiles
    currentfile = sprintf('setB%d.mat',k);
    Dfile = load(['C:\Users\isabe\Documents\MATLAB\Master\Data\Box_data\',...
        currentfile]);
    dataset{k}=Dfile.table;
end

%Using the function FFTset.m on all twenty data sets, storing the data for
%each separately:
TF = cell(1,20);
DF_nom = cell(1,20);
DF_denom = cell(1,20);
DI = cell(1,20);
ka = cell(1,20);
fvec_ka = cell(1,20);

%Loading file containing calculated sound velocity:
Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\Data\Box_data\sound_vel_data_box.mat');
cdata = Dfile.sound_vel_data;

cdata_avg = zeros(20,1);
%Calculating average velocity of sound for each measurement series:
for k = 1:20
     cdata_avg(k,:) = mean(cdata(k,:),2); %Temperature
end

%Using box length (average for all sides):
radius = 349.7341667/1000;

%Perfoming Fourier Transformation on Impulse Responses:
for k = 1:20
    if k<=4
        lspnum = 2;
    else
        lspnum = 1;
    end
    [TF{k},DF_nom{k},DF_denom{k},DI{k},ka{k},fvec_ka{k}]= ...
        FFTset_v4_diffwindow(dataset{k},lspnum,cdata_avg(k,:),radius);
end

%% Background noise 

%Loading data file for background noise (ventillation off):
bg_noise_1 = load('C:\Users\isabe\Documents\MATLAB\Master\Data\Box_data\Box_noise_Thur.mat');
bg_noise_1 = bg_noise_1.table;
data = zeros(20000,2);
data(:,1) = bg_noise_1;

%Noise for ventillation on:
bg_noise_2 = load('C:\Users\isabe\Documents\MATLAB\Master\Data\Box_data\Box_noise_vent_1.mat');
bg_noise_2 = bg_noise_2.table;
data(:,2) = bg_noise_2;

%Defining start and end point of the data to be considered in the FFT:
w_start = 250;
w_end = 550;

%Calculating FFT with a window:
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

TF_bg = fft(data(w_start:w_end,:).*full_win,nfft);

%Discarding the mirrored negative frequencies:
TF_bg = TF_bg(1:nfft/2,:);

%% Figures

% %Figure containing the directivity index for all eight measurement series:
% figure(1)
% semilogx(ka_1,DI_1,'-',ka_2,DI_2,'-',ka_3,DI_3,'-',ka_4,DI_4,'-',...
%     ka_5,DI_5,'--',ka_6,DI_6,'--',ka_7,DI_7,'--',ka_8,DI_8,'--')
% xlim([0.1 60])
% g = legend('1st','2nd','3rd','4th','5th','6th','7th','8th');
% set(g,'Location','best')
% title({' ','Directivity index of measurement','sets 1 through 8',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude',' '})
% set(gcf,'color','w')
% grid on

% %Figure plotting the difference in directivity index for the measurement
% %series done for loudspeaker #1:
% figure(2)
% semilogx(ka_1,DI_diff_1,'-o',ka_1,DI_diff_2,'-o')
% xlim([0.5 60])
% g = legend('Diff. sets 1 & 2','Diff. sets 3 & 4');
% set(g,'Location','best')
% title({' ','Difference in DI for,','loudspeaker #1',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude',' '})
% set(gcf,'color','w')
% grid on


%Figure containing the numerator and denominator plotted separately for
%each measurement series for loudspeaker #1, with the denominator shifted
%upward by a factor for comparability:

% %Factor for shifting denominator up to the magnitude of the numerator:
% f1 = 1;%1.98;
% figure(3)
% loglog(ka_1,DF_nom_1,'-',ka_1,DF_denom_1*f1,'--',ka_2,DF_nom_2,'-',...
%     ka_2,DF_denom_2*f1,'--',ka_3,DF_nom_3,'-',ka_3,DF_denom_3*f1,'--',...
%     ka_4,DF_nom_4,'-',ka_4,DF_denom_4*f1,'--',ka_1,abs(TF_bg(:,1)),'-.',ka_1,abs(TF_bg(:,2)),'-.')
% g = legend('1st num','1st denom','2nd num','2nd denom','3rd num',...
%     '3rd denom','4th num','4th denom','noise floor vent. off','noise floor vent. on');
% set(g,'Location','best')
% title({' ','Plot of numerators and denominators, sets',...
%     '1 through 4 with the background noise level',...
%     '(denominator shifted up)',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude',' '})
% set(gcf,'color','w')
% grid on


% %Figure containing the numerator and denominator plotted separately for
% %each measurement series for loudspeaker #2, with the denominator shifted
% %upward by a factor for comparability:
% %Factor for shifting denominator up to the magnitude of the numerator:
% f2 = 1;%2.13;
% figure(4)
% loglog(ka_5,DF_nom_5,'-',ka_5,DF_denom_5*f2,'--',ka_6,DF_nom_6,'-',...
%     ka_6,DF_denom_6*f2,'--',ka_7,DF_nom_7,'-',ka_7,DF_denom_7*f2,'--',...
%     ka_8,DF_nom_8,'-',ka_8,DF_denom_8*f2,'--')
% g = legend('5th num','5th denom','6th num','6th denom','7th num',...
%     '7th denom','8th num','8th denom');
% set(g,'Location','best')
% title({' ','Plot of numerators and denominators,',...
%     'sets 5 through 8 (denom. shifted up)',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude',' '})
% set(gcf,'color','w')
% grid on

% %Signal to noise ratio for both loudspeakers (from measurement series 1&5):
% figure(5)
% semilogx(ka_1,10*log10(DF_nom_1(:)./abs(TF_bg(:,1).^2)),'-',...
%     ka_5,10*log10(DF_nom_5(:)./abs(TF_bg(:,1).^2)),'-.')
% title({' ','Signal to noise ratio for',...
%     'loudspeaker #1 & #2', '(data from series 1 & 5)',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude, dB',' '})
% set(gcf,'color','w')
% g = legend('LS #2','LS #1');
% set(g,'Location','best')
% grid on

% %Signal to noise ratio for loudspeaker #1 (from measurement series 1):
% figure(6)
% semilogx(ka_1,10*log10(DF_nom_1(:)./abs(TF_bg(:).^2)),'-')
% title({' ','Signal to noise ratio for','loudspeaker #1 (series 1)',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude, dB',' '})
% set(gcf,'color','w')
% grid on
% 
% %Signal to noise ratio for loudspeaker #2 (from measurement series 5):
% figure(7)
% semilogx(ka_5,10*log10(DF_nom_5(:)./abs(TF_bg(:).^2)),'-')
% title({' ','Signal to noise ratio for','loudspeaker #2 (series 5)',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude, dB',' '})
% set(gcf,'color','w')
% grid on

% %Plot of numerator and denominator for one measurement series with each
% %loudspeaker, as well as noise floor w/vent. on & off
% figure(6)
% fig6 = loglog(ka{5},DF_nom{5},'-',ka{5},DF_denom{5},'--',...
%     ka{1},DF_nom{1},'-',ka{1},DF_denom{1},'--',...
%     ka{18},abs(TF_bg(:,1)),'-.',ka{1},abs(TF_bg(:,2)),'-.'); %Is the choice in ka correct here?
% set(fig6,'Linewidth',1)
% fig6(1).Color = [0 0.4470 0.7410]; %blue
% fig6(2).Color = [0 0.4470 0.7410];
% fig6(3).Color = [0.8500 0.3250 0.0980]; %red
% fig6(4).Color = [0.8500 0.3250 0.0980];
% fig6(5).Color = [0.9290 0.6940 0.1250]; %yellow
% fig6(6).Color = [0.4660 0.6740 0.1880]; %green
% g = legend('LS#1 numerator','LS#1 denominator','LS#2 numerator',...
%     'LS#2 denominator','noise floor vent. off','noise floor vent. on');
% set(g,'Location','best')
% % title({'Plot of numerators and denominators, sets',...
% %     '1 and 5 with the background noise levels',''})
% xlabel({'kl, [1]'})
% ylabel({'Magnitude, [1]'})
% set(gcf,'color','w')
% grid on

AA = ((max([max(10*log10(abs(TF_bg(35:end,2)))) (max(10*log10(abs(TF_bg(35:end,1)))))])) + 30);
BB = AA - 30;
SNR_y = AA*ones(1,100);
SNR_x = linspace(ka{1}(1,2),ka{1}(1,end),100);
Noisefloor_y = BB*ones(1,100);

figure(7)
fig7 = semilogx(ka{5},10*log10(DF_nom{5}),'-',ka{5},...
    10*log10(DF_denom{5}),'-.',ka{1},10*log10(DF_nom{1}),'-',...
    ka{1},10*log10(DF_denom{1}),'-.',ka{18},10*log10(abs(TF_bg(:,1))),...
    '-.',ka{1},10*log10(abs(TF_bg(:,2))),'-.',SNR_x,Noisefloor_y,'--',...
    SNR_x,SNR_y,'-'); %Is the choice in ka correct here?
set(fig7,'Linewidth',1)
fig7(1).Color = [0 0.4470 0.7410]; %blue
fig7(2).Color = [0 0.4470 0.7410];
fig7(3).Color = [0.8500 0.3250 0.0980]; %red
fig7(4).Color = [0.8500 0.3250 0.0980];
fig7(5).Color = [0.9290 0.6940 0.1250]; %yellow
fig7(6).Color = [0.4940 0.1840 0.5560]; %purple
fig7(8).Color = [0.4660 0.6740 0.1880]; %green
fig7(7).Color = [0.9290 0.6940 0.1250]; %yellow
xlim([0.03 165])
ylim([-30 70])
g = legend('LS#1 numerator','LS#1 denominator','LS#2 numerator',...
    'LS#2 denominator','Noise vent. off','Noise vent. on',...
    'Noise floor','Noise floor +30 dB');
set(g,'Location','best')
% title({'Plot of numerators and denominators, sets',...
%     '1 and 5 with the background noise levels',''})
xlabel({'kl, [1]'})
ylabel({'Magnitude, [dB]'})
set(gcf,'color','w')
grid on

%% Lower cutoff values from SNR = 30 dB:
%LS1: 1.699823433909710 ca. 267 Hz
%LS2: 0.747671418147417 ca. 117 Hz
%LS2: 139.4220276990395 upper limit ca. 21852.53906250000 Hz

