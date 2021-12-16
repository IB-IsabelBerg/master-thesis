% Script to plot the impulse response for measurements done on the cylinder
% and on the box, as well as plot the window function used when performing
% the discrete Fourier transform on the impulse response during processing.
%
% Written 29.10.21 by Isabel Berg

clear all
close all

%% Load data files

cylfile = load('C:\Users\isabe\Documents\MATLAB\Master\Data\data_set_5.mat');
boxfile = load('C:\Users\isabe\Documents\MATLAB\Master\Data\Box_data\setB1.mat');
cyldata = cylfile.table;
boxdata = boxfile.table;
boxdata = (-1)*boxdata; %Measurement series was done with poles switched

%% Do windowing
%Without edge reflection:
w_start = 250;
w_end = 600;

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

% %Making half Hanning for start as well:
% n_start = 10;
% start_win = hanning(n_start*2);
% start_win = start_win(1:(n_start));
% full_win(1:(n_start)) = start_win(:);

x_win = (w_start:w_end);

max_cyldata = (max(max(cyldata)));
max_boxdata = (max(max(boxdata)));

for k = 1:size(cyldata,2)
    normalized_cyldata(:,k) = cyldata(:,k)./max_cyldata;
    normalized_boxdata(:,k) = boxdata(:,k)./max_boxdata;
end

x_data = linspace(1,20000,20000);

% TF_set_cyl = fft(cyldata(w_start:w_end,:).*full_win,nfft);
% TF_set_box = fft(boxdata(w_start:w_end,:).*full_win,nfft);
% 
% %Discarding mirrored half:
% TF_set_cyl = TF_set_cyl(1:nfft/2,:);
% TF_set_box = TF_set_box(1:nfft/2,:);

%% Figures

%Plot of IR
% figure(1)
% plot(cyldata)
% grid on
% % 
% figure(2)
% plot(full_win)
% grid on

% figure(3)
% plot(x_data,normalized_cyldata(:,91),'-',x_win,full_win,'-',x_data,normalized_cyldata(:,1),'--');
% grid on
% xlim([0 700]);
% ylim([-0.1 1.1]);
% 
% figure(4)
% plot(normalized_cyldata)
% grid on
x = (1:1:length(normalized_cyldata));

figure(5);
fig5 = plot(x,normalized_cyldata,x_win,full_win);
set(fig5,'LineWidth',1)
grid on
%hold on
%title({'Cylinder data with half Hanning window',''})
ylabel({'Impulse response (normalized), [1]'})
xlabel({'Sample no., [1]'})
%fig1 = plot(x_win,full_win);
%set(fig1,'LineWidth',1)
g = legend('Half Hanning window');
set(g,'Location','best');
xlim([0 700]);
ylim([-0.75 1.1]);
set(gcf,'color','w');
hold off

figure(6)
fig6 = plot(x,normalized_boxdata,x_win,full_win);
set(fig6,'LineWidth',1)
grid on
%hold on
%title({'Box data with half Hanning window',''})
ylabel({'Impulse response (normalized), [1]'})
xlabel({'Sample no., [1]'})
%fig2 = plot(
%set(fig2,'LineWidth',1)
g = legend('Half Hanning window');
set(g,'Location','best');
xlim([0 700]);
ylim([-0.75 1.1]);
set(gcf,'color','w');
hold off

