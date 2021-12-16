% Script that plots the measured directivity index for the cylinder, along
% with a plot of the difference between the directivity index for the two
% different measurement setups. Compares both the average of the DI with
% maximum deviation, as well as the difference between the max deviations.

% Written 28.11.21 by Isabel Berg
% Modified 01.12.21 to use average ka values


%%%, and best index??
clear all
close all


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


%% DF and DI

%Calculating Directivity Factor:
DF = cell(1,numfiles);
for k = 1:numfiles
    DF{k} = DF_nom{k}./DF_denom{k};
end

DF_avg = cell(1,2);
ka_avg = cell(1,2);
%Averaging over four series for each setup:
for l = 1:2
    if l == 1
        DF_avg{l} = (DF{1} + DF{2} + DF{3} + DF{4})/4;
        ka_avg{l} = (ka{1} + ka{2} + ka{3} + ka{4})/4;
    else
        DF_avg{l} = (DF{5} + DF{6} + DF{7} + DF{8})/4;
        ka_avg{l} = (ka{5} + ka{6} + ka{7} + ka{8})/4;
    end
end

%% Added 01.12
%Best index for averaging:
indexnumber_in_ka = cell(1,8);
for l = 1:4
    indexnumber_in_ka{l} = findbestkavecindices(ka{l},ka_avg{1});
    indexnumber_in_ka{4+l} = findbestkavecindices(ka{4+l},ka_avg{2});
end
% 
DF_avg_best_index = cell(1,2);
for m = 1:length(DF{1})
    for l = 1:2
        if l == 1
            DF_avg_best_index{l}(m) = (DF{1}(indexnumber_in_ka{1}(m))...
                + DF{2}(indexnumber_in_ka{2}(m))...
                + DF{3}(indexnumber_in_ka{3}(m))...
                + DF{4}(indexnumber_in_ka{4}(m)))/4;
        else
            DF_avg_best_index{l}(m) = (DF{5}(indexnumber_in_ka{5}(m))...
                + DF{6}(indexnumber_in_ka{6}(m))...
                + DF{7}(indexnumber_in_ka{7}(m))...
                + DF{8}(indexnumber_in_ka{8}(m)))/4;
        end
    end
end

diff_DF_avg = cell(1,2);
diff_DF_avg{1} = DF_avg_best_index{1} - DF_avg{1};
diff_DF_avg{2} = DF_avg_best_index{2} - DF_avg{2};

%%

DF_dev = cell(1,2);
DF_max_dev = cell(1,2);
DF_min_dev = cell(1,2);
%Calculating max and min deviations from average within each set:
for l = 1:2
    if l == 1
        DF_dev{l} = [(DF{1}); DF{2}; DF{3}; DF{4}];
    else
        DF_dev{l} = [(DF{5}); DF{6}; DF{7}; DF{8}];
    end
    DF_max_dev{l} = max(DF_dev{l},[],1);
    DF_min_dev{l} = min(DF_dev{l},[],1);
end

DI_avg = cell(1,2);
DI_upper_dev = cell(1,2);
DI_lower_dev = cell(1,2);
%Turning Directivity Factor into Directivity Index:
for l = 1:2
    DI_avg{l} = 10*log10(DF_avg{l});
    DI_upper_dev{l} = 10*log10(DF_max_dev{l});
    DI_lower_dev{l} = 10*log10(DF_min_dev{l});
end

%% Comparison of DI; difference

%S2 minus S1
DI_mean_diff_A = DI_avg{2} - DI_avg{1};
DI_max_diff_A = DI_upper_dev{2} - DI_lower_dev{1};

ka_avg_diff = (ka_avg{1} + ka_avg{2})/2;

%% Best index
DF_dev_best_index = cell(1,2);
DF_max_dev_best_index = cell(1,2);
DF_min_dev_best_index = cell(1,2);
%Calculating max and min deviations from average within each set:
for m = 1:length(DF{1})
    for l = 1:2
        if l == 1
            DF_dev_best_index{l}(1:4,m) = [DF{1}(indexnumber_in_ka{1}(m)); ...
                DF{2}(indexnumber_in_ka{2}(m)); ...
                DF{3}(indexnumber_in_ka{3}(m)); ...
                DF{4}(indexnumber_in_ka{4}(m))];
        elseif l == 2
            DF_dev_best_index{l}(1:4,m) = [DF{5}(indexnumber_in_ka{5}(m)); ...
                DF{6}(indexnumber_in_ka{6}(m)); ...
                DF{7}(indexnumber_in_ka{7}(m)); ...
                DF{8}(indexnumber_in_ka{8}(m))];
        end
        DF_max_dev_best_index{l}(m) = max(DF_dev_best_index{l}(:,m),[],1);
        DF_min_dev_best_index{l}(m) = min(DF_dev_best_index{l}(:,m),[],1);
    end
end

DI_avg_best_index = cell(1,2);
DI_upper_dev_best_index = cell(1,2);
DI_lower_dev_best_index = cell(1,2);
%Turning Directivity Factor into Directivity Index:
for l = 1:2
    DI_avg_best_index{l} = 10*log10(DF_avg_best_index{l});
    DI_upper_dev_best_index{l} = 10*log10(DF_max_dev_best_index{l});
    DI_lower_dev_best_index{l} = 10*log10(DF_min_dev_best_index{l});
end

% Comparison of DI; difference

%S1 and S2
DI_mean_diff_A_best_index = DI_avg_best_index{2} - DI_avg_best_index{1};
DI_max_diff_A_best_index = DI_upper_dev_best_index{2} - ...
    DI_lower_dev_best_index{1};

diff_DI_avg = cell(1,2);
diff_DI_avg{1} = DI_avg_best_index{1} - DI_avg{1};
diff_DI_avg{2} = DI_avg_best_index{2} - DI_avg{2};

%% Figures

%Difference between same index and best index
figure(3)
semilogx(ka_avg_diff,(DI_mean_diff_A-DI_mean_diff_A_best_index),...
    ka_avg_diff,(DI_max_diff_A-DI_max_diff_A_best_index))
xlim([0.435441101586095 43.853587469566044])
xlabel({'ka, [1]',' '})
g = legend('Difference, mean','Difference, max');
ylabel({'Difference in directivity index, [dB]'})
set(gcf,'color','w');
grid on
title('Difference between same index and best index')

% Directivity index
figure(7)
fig6 = semilogx(ka_avg{1},DI_avg_best_index{1},'-',...
    ka_avg{2},DI_avg_best_index{2},'-',...
    ka_avg{1},DI_upper_dev_best_index{1},'-.',...
    ka_avg{1},DI_lower_dev_best_index{1},'-.',...
    ka_avg{2},DI_upper_dev_best_index{2},'-.',...
    ka_avg{2},DI_lower_dev_best_index{2},'-.');
set(fig6,'Linewidth',1)
fig6(1).Color = [0 0.4470 0.7410]; %blue
fig6(2).Color = [0.8500 0.3250 0.0980]; %red
fig6(3).Color = [0 0.4470 0.7410]; %blue
fig6(4).Color = [0 0.4470 0.7410]; %blue
fig6(5).Color = [0.8500 0.3250 0.0980]; %red
fig6(6).Color = [0.8500 0.3250 0.0980]; %red
xlim([0.435441101586095 43.853587469566044])
ylim([1 8]);
xlabel({'ka, [1]'})
ylabel({'Directivity index, [dB]'})
set(gcf,'color','w');
g = legend(...
    'LS#1',...
    'LS#2');
set(g,'Location','best');
grid on

% Difference in directivity index
figure(15)
fig15 = semilogx(ka_avg_diff,DI_mean_diff_A_best_index,'-',...
    ka_avg_diff,DI_max_diff_A_best_index,'-.');
set(fig15,'Linewidth',1)
xlim([0.435441101586095 43.853587469566044])
ylim([-0.3 1.5])
fig15(1).Color = [0 0.4470 0.7410]; %blue
fig15(2).Color = [0 0.4470 0.7410];
xlabel({'ka, [1]',' '})
ylabel({'Difference in directivity index, [dB]'})
set(gcf,'color','w');
g = legend('LS#2 vs. LS#1 (mean)','LS#2 vs. LS#1 (max)');
grid on

% ka value problem illustration
sn = 4561;
figure(1)
fig1 = semilogx(ka_avg{1},DF{1},'-o',ka_avg{1},DF{2},'-o',...
    ka_avg{1},DF{3},'-o',ka_avg{1},DF{4},'-o',ka_avg{1},DF_avg{1},'-o',...
    [ka_avg{1}(sn) ka_avg{1}(sn) ka_avg{1}(sn) ka_avg{1}(sn) ka_avg{1}(sn)],...
    [DF{1}(sn) DF{2}(sn) DF{3}(sn) DF{4}(sn) DF_avg{1}(sn)],...
    '*k');
set(fig1,'Linewidth',1)
xlim([24.5 24.57])
ylim([2.69 2.77])
xlabel({'ka, [1]'})
ylabel({'Directivity factor, [1]'})
% ylim([2.755 2.795])
% xlim([24.77 24.83])
% ylim([2.755 2.795])
grid on
g = legend('Series 1','Series 2','Series 3','Series 4','Average');
set(g,'Location','best')
set(gcf,'color','w');

figure(2)
fig2 = semilogx(ka{1},DF{1},'-o',ka{2},DF{2},'-o',ka{3},DF{3},'-o',...
    ka{4},DF{4},'-o',ka_avg{1},DF_avg{1},'-o',...
    [ka{1}(sn) ka{2}(sn) ka{3}(sn) ka{4}(sn) ka_avg{1}(sn)],...
    [DF{1}(sn) DF{2}(sn) DF{3}(sn) DF{4}(sn) DF_avg{1}(sn)],...
    '*k');
set(fig2,'Linewidth',1)
xlim([24.5 24.57])
ylim([2.69 2.77])
xlabel({'ka, [1]'})
ylabel({'Directivity factor, [1]'})
% xlim([24.77 24.83])
% ylim([2.755 2.795])
grid on
g = legend('Series 1','Series 2','Series 3','Series 4','Average');
set(g,'Location','best')
set(gcf,'color','w');


% Test, diff for best index
% figure(4)
% semilogx(ka_avg{1},diff_DI_avg{1},ka_avg{2},diff_DI_avg{2})
% grid on
% xlim([0.435441101586095 43.853587469566044])
% xlabel({'ka, [1]',' '})
% ylabel({'Difference in directivity index, [dB]'})
% set(gcf,'color','w');

% Before best index
% figure(14)
% fig141 = semilogx(ka_avg_diff,DI_mean_diff_A,'-',ka_avg_diff,DI_max_diff_A,...
%     '-.');
% set(fig141,'Linewidth',1)
% xlim([0.435441101586095 43.853587469566044])
% ylim([-0.3 1.5])
% fig141(1).Color = [0 0.4470 0.7410]; %blue
% fig141(2).Color = [0 0.4470 0.7410];
% xlabel({'ka, [1]',' '})
% ylabel({'Difference in directivity index, [dB]'})
% set(gcf,'color','w');
% g = legend('LS#1 vs. LS#2 (mean)','LS#1 vs. LS#2 (max)');
% grid on

%% Lower cutoff from SNR = 30 dB:
%LS1: 0.435441101586095 ca. 237 Hz
%LS2: 0.204445629228746 ca. 111 Hz
%LS2: 43.853587469566044 ca. 23879.88281250000 Hz upper limit



