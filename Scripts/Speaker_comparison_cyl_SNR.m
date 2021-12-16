%Creates plot of the directivity factor, but with the numerator and the
%denominator plotted separately. Includes the noise floor. Utilizes the
%function FFTset_v4.m, as well as importing measurement, noise, and C data
%files. Written for the cylinder scattering object.

%Original written 12.05.21 by Isabel Berg
%Updated 19.10.21 to use FFTset_v4 (including C data and radius), and make
%the code more efficient.
clear all
close all
%% REQUIRES THE FUNCTION "FFTset_v4.m"

%% REQUIRES LOCALLY STORED MEASUREMENT, NOISE, AND C DATA FILES (TOT: 10)

%% Opening data sets, using function for running FFT and calculating DF
dataset = cell(1,8);
%Loading all data sets (PLACEMENTS TO BE CHANGED WHEN RUN ELSEWHERE):
numfiles = 8;
for k = 1:numfiles
    currentfile = sprintf('data_set_%d.mat',k);
    Dfile = load(['C:\Users\isabe\Documents\MATLAB\Master\Data\',...
        currentfile]);
    dataset{k}=Dfile.table;
end

%Using the function FFTset.m on all twenty data sets, storing the data for
%each separately:
TF = cell(1,8);
DF_nom = cell(1,8);
DF_denom = cell(1,8);
DI = cell(1,8);
ka = cell(1,8);
fvec_ka = cell(1,8);

%Loading file containing calculated sound velocity:
Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\Data\C_data_cyl.mat');
cdata = Dfile.C_data;

cdata_avg = zeros(8,1);
%Calculating average velocity of sound for each measurement series:
for k = 1:8
     cdata_avg(k,:) = mean(cdata(k,:),2); %Temperature
end

%Cylinder radius (measured):
radius = (0.20043)/2; % [m]

%Perfoming Fourier Transformation on Impulse Responses:
for k = 1:8
    if k>=4
        lspnum = 2;
    else
        lspnum = 1;
    end
    [TF{k},DF_nom{k},DF_denom{k},DI{k},ka{k},fvec_ka{k}]= ...
        FFTset_v4(dataset{k},lspnum,cdata_avg(k,:),radius);
end


%% Background noise 

%Loading data file for background noise:
bg_noise = load('C:\Users\isabe\Documents\MATLAB\Master\Data\bg_noise.mat');
bg_noise = bg_noise.table;
data = bg_noise;

%Defining start and end point of the data to be considered in the FFT:
w_start = 250;
w_end = 600;

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
% 
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
% 
% 
% %Figure containing the numerator and denominator plotted separately for
% %each measurement series for loudspeaker #1, with the denominator shifted
% %upward by a factor for comparability:
% 
% %Factor for shifting denominator up to the magnitude of the numerator:
% f1 = 1;%1.98;
% figure(3)
% loglog(ka_1,DF_nom_1,'-',ka_1,DF_denom_1*f1,'--',ka_2,DF_nom_2,'-',...
%     ka_2,DF_denom_2*f1,'--',ka_3,DF_nom_3,'-',ka_3,DF_denom_3*f1,'--',...
%     ka_4,DF_nom_4,'-',ka_4,DF_denom_4*f1,'--',ka_1,abs(TF_bg),'-.')
% g = legend('1st num','1st denom','2nd num','2nd denom','3rd num',...
%     '3rd denom','4th num','4th denom','noise floor');
% set(g,'Location','best')
% title({' ','Plot of numerators and denominators, sets',...
%     '1 through 4 with the background noise level',...
%     '(denominator shifted up)',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude',' '})
% set(gcf,'color','w')
% grid on
% 
% 
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
% semilogx(ka_1,10*log10(DF_nom_1(:)./abs(TF_bg(:).^2)),'-',...
%     ka_5,10*log10(DF_nom_5(:)./abs(TF_bg(:).^2)),'-.')
% title({' ','Signal to noise ratio for',...
%     'loudspeaker #1 & #2', '(data from series 1 & 5)',''})
% xlabel({'ka',' '})
% ylabel({' ','Magnitude, dB',' '})
% set(gcf,'color','w')
% g = legend('LS #1','LS #2');
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
% figure(7)
% fig7 = loglog(ka{1},DF_nom{1},'-',ka{1},DF_denom{1},'--',ka{5},...
%     DF_nom{5},'-',ka{5},DF_denom{5},'--',...
%     ka{8},abs(TF_bg(:,1)),'-.'); %Is the choice in ka correct here? Check when bg noise measured
% set(fig7,'Linewidth',1)
% fig7(1).Color = [0 0.4470 0.7410]; %blue
% fig7(2).Color = [0 0.4470 0.7410];
% fig7(3).Color = [0.8500 0.3250 0.0980]; %red
% fig7(4).Color = [0.8500 0.3250 0.0980];
% fig7(5).Color = [0.9290 0.6940 0.1250]; %yellow
% g = legend('LS#1 numerator','LS#1 denominator','LS#2 numerator',...
%     'LS#2 denominator','noise floor vent. off');
% set(g,'Location','best')
% % title({'Plot of numerators and denominators, sets',...
% %     '1 and 5 with the background noise levels',''})
% xlabel({'ka, [1]'})
% ylabel({'Magnitude, [1]'})
% set(gcf,'color','w')
% grid on
AA = (max(10*log10(abs(TF_bg(:,1)))) + 30);
BB =  AA-30;
SNR_y = AA*ones(1,100);
SNR_x = linspace(ka{1}(1,2),ka{1}(1,end),100);
Noisefloor_y = BB*ones(1,100);

figure(8)
fig8 = semilogx(ka{1},10*log10(DF_nom{1}),'-',ka{1},10*log10(DF_denom{1}),'-.',ka{5},...
    10*log10(DF_nom{5}),'-',ka{5},10*log10(DF_denom{5}),'-.',...
    ka{8},10*log10(abs(TF_bg(:,1))),'-.',SNR_x,Noisefloor_y,'--',...
    SNR_x,SNR_y,'-'); %Is the choice in ka correct here? Check when bg noise measured
set(fig8,'Linewidth',1)
xlim([0.01 50])
ylim([-30 70])
fig8(1).Color = [0 0.4470 0.7410]; %blue
fig8(2).Color = [0 0.4470 0.7410];
fig8(3).Color = [0.8500 0.3250 0.0980]; %red
fig8(4).Color = [0.8500 0.3250 0.0980];
fig8(5).Color = [0.9290 0.6940 0.1250]; %yellow
fig8(7).Color = [0.4660 0.6740 0.1880]; %green
fig8(6).Color = [0.9290 0.6940 0.1250]; %green
g = legend('LS#1 numerator','LS#1 denominator','LS#2 numerator',...
    'LS#2 denominator','Noise vent. off','Noise floor',...
    'Noise floor +30 dB');
set(g,'Location','best')
% title({'Plot of numerators and denominators, sets',...
%     '1 and 5 with the background noise levels',''})
xlabel({'ka, [1]'})
ylabel({'Magnitude, [dB]'})
set(gcf,'color','w')
grid on

%% Lower cutoff from SNR = 30 dB:
%LS1: 0.435441101586095 ca. 237 Hz
%LS2: 0.204445629228746 ca. 111 Hz
%LS2: 43.853587469566044 ca. 23879.88281250000 Hz upper limit

% figure(9)
% grid on
% T = ka{1};
% S = fvec;
% D = 10*log10(DF_nom{1});
% x1 = T;
% y1 = D;
% x2 = S;
% y2 = D;
% xlabels{1} = 'Temperature (C)';
% xlabels{2} = 'Salinity';
% ylabels{1} = 'Depth(m)';
% ylabels{2} = 'Depth(m)';
% hl1=plot(x1,y1,'Color','w');
% ax(1)=gca;
% set(ax(1),'Position',[0.12 0.12 0.75 0.70])
% set(ax(1),'XColor','k','YColor','k','XScale','log');
% ax(2)=axes('Position',get(ax(1),'Position'),...
%    'XAxisLocation','top',...
%    'YAxisLocation','right',...
%    'Color','none',...
%    'XColor','k','YColor','w','XScale','log');
% set(ax(2),'XScale','log');
% set(ax,'box','off')
% hl2=line(x2,y2,'Color','k','Parent',ax(2));
% 
% %label the two x-axes
% set(get(ax(1),'xlabel'),'string',xlabels{1})
% set(get(ax(2),'xlabel'),'string',xlabels{2})
% set(get(ax(1),'ylabel'),'string',ylabels{1})
% set(get(ax(2),'ylabel'),'string',ylabels{2})


