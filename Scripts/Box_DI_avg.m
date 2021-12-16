%Script that calculates average Directivity Index for measurements
%performed on the box object, including max-min deviation from calculated
%average.

%Written 16.09.21 by Isabel Berg
%Edited 20.09.21 for implementation of sound velocity data calculated from
%monitored variables (temperature, humidity, pressure) during measurements
%Edited 26.10.21 for plots of difference in DI
%Edited 06.12.21 to use bestindex

clear all

%% To do

% Clean up plots; what needs to be included in the results?

%Compare plots with separate ka vectors (displacement) with those of
%average ka vec. Figure 2 -- set xlim and ylim for better comparison

%Create new plots, with difference in DI plotted. 

%Averaging of fvec in the event of plotting the averages of the different
%measurement setups with f axis instead of ka axis? Or just settle for
%using a "kl" axis?

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
fvec = cell(1,20);


%Loading file containing calculated sound velocity:
Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\Data\Box_data\sound_vel_data_box.mat');
cdata = Dfile.sound_vel_data;

cdata_avg = zeros(20,1);
%Calculating average velocity of sound for each measurement series:
for k = 1:20
     cdata_avg(k,:) = mean(cdata(k,:),2); %Temperature
end

radius = 349.7341667/1000;

%Perfoming Fourier Transformation on Impulse Responses:
for k = 1:20
    if k<=4
        lspnum = 2;
    else
        lspnum = 1;
    end
    [TF{k},DF_nom{k},DF_denom{k},DI{k},ka{k},fvec{k}]= ... % technically kl
        FFTset_v4(dataset{k},lspnum,cdata_avg(k,:),radius);
end

%% DF and DI

%Calculating Directivity Factor:
DF = cell(1,20);
for k = 1:20
    DF{k} = DF_nom{k}./DF_denom{k};
end

DF_avg = cell(1,5);
ka_avg = cell(1,5);
%Averaging over four series for each setup:
for l = 1:5
    if l == 2
        DF_avg{l} = (DF{5} + DF{10} + DF{11} + DF{12})/4;
        ka_avg{l} = (ka{5} + ka{10} + ka{11} + ka{12})/4;
    elseif l == 3
        DF_avg{l} = (DF{6} + DF{7} + DF{8} + DF{9})/4;
        ka_avg{l} = (ka{6} + ka{7} + ka{8} + ka{9})/4;
    else
        DF_avg{l} = (DF{4*(l-1)+1} + DF{4*(l-1)+2} + DF{4*(l-1)+3} + ...
            DF{4*(l-1)+4})/4;
        ka_avg{l} = (ka{4*(l-1)+1} + ka{4*(l-1)+2} + ka{4*(l-1)+3} + ...
            ka{4*(l-1)+4})/4;
    end
end

%% Added 06.12
%Best index for averaging:
indexnumber_in_ka = cell(1,20);
for l = 1:5
    if l == 2
        indexnumber_in_ka{5} = findbestkavecindices(ka{5},ka_avg{1});
        indexnumber_in_ka{10} = findbestkavecindices(ka{10},ka_avg{1});
        indexnumber_in_ka{11} = findbestkavecindices(ka{11},ka_avg{1});
        indexnumber_in_ka{12} = findbestkavecindices(ka{12},ka_avg{1});
    elseif l == 3
        for m = 1:4
            indexnumber_in_ka{(5+m)} = ...
                findbestkavecindices(ka{(5+m)},ka_avg{1});
        end
    else
        for m = 1:4
            indexnumber_in_ka{(4*(l-1)+m)} = ...
                findbestkavecindices(ka{(4*(l-1)+m)},ka_avg{1});
        end
    end
end

DF_avg_best_index = cell(1,5);
for m = 1:length(DF{1})
    for l = 1:5
        if l == 2
            DF_avg_best_index{l}(m) = (DF{5}(indexnumber_in_ka{5}(m))...
                + DF{10}(indexnumber_in_ka{10}(m))...
                + DF{11}(indexnumber_in_ka{11}(m))...
                + DF{12}(indexnumber_in_ka{12}(m)))/4;
        elseif l == 3
            DF_avg_best_index{l}(m) = (DF{6}(indexnumber_in_ka{6}(m))...
                + DF{7}(indexnumber_in_ka{7}(m))...
                + DF{8}(indexnumber_in_ka{8}(m))...
                + DF{9}(indexnumber_in_ka{9}(m)))/4;
        else
            DF_avg_best_index{l}(m) = ...
                (DF{(4*(l-1)+1)}(indexnumber_in_ka{(4*(l-1)+1)}(m))...
                + DF{(4*(l-1)+2)}(indexnumber_in_ka{(4*(l-1)+2)}(m))...
                + DF{(4*(l-1)+3)}(indexnumber_in_ka{(4*(l-1)+3)}(m))...
                + DF{(4*(l-1)+4)}(indexnumber_in_ka{(4*(l-1)+4)}(m)))/4;
        end
    end
end

%%
DF_dev = cell(1,5);
DF_max_dev = cell(1,5);
DF_min_dev = cell(1,5);
%Calculating max and min deviations from average within each set:
for l = 1:5
    if l == 2
        DF_dev{l} = [(DF{5}); DF{10}; DF{11}; DF{12}];
    elseif l == 3
        DF_dev{l} = [(DF{6}); DF{7}; DF{8}; DF{9}];
    else
        DF_dev{l} = [(DF{4*(l-1)+1}); DF{4*(l-1)+2}; DF{4*(l-1)+3}; ...
            (DF{4*(l-1)+4})]; 
    end
    DF_max_dev{l} = max(DF_dev{l},[],1);
    DF_min_dev{l} = min(DF_dev{l},[],1);
end

DI_avg = cell(1,5);
DI_upper_dev = cell(1,5);
DI_lower_dev = cell(1,5);
%Turning Directivity Factor into Directivity Index:
for l = 1:5
    DI_avg{l} = 10*log10(DF_avg{l});
    DI_upper_dev{l} = 10*log10(DF_max_dev{l});
    DI_lower_dev{l} = 10*log10(DF_min_dev{l});
end


%% Comparison of DI; difference

%S1 and S2
DI_mean_diff_A = DI_avg{1} - DI_avg{2};
DI_max_diff_A = DI_upper_dev{1} - DI_lower_dev{2};
ka_avg_A = (ka_avg{2} + ka_avg{1})/2;

%Ventilation
DI_mean_diff_B = DI_avg{2} - DI_avg{3};
DI_max_diff_B = DI_upper_dev{2} - DI_lower_dev{3};
ka_avg_B = (ka_avg{2} + ka_avg{3})/2;

%Adjusted mic
DI_mean_diff_C = DI_avg{2} - DI_avg{4};
DI_max_diff_C = DI_upper_dev{2} - DI_lower_dev{4};
ka_avg_C = (ka_avg{2} + ka_avg{4})/2;

%Adjusted loudspeaker
DI_mean_diff_D = DI_avg{2} - DI_avg{5};
DI_max_diff_D = DI_upper_dev{2} - DI_lower_dev{5};
ka_avg_D = (ka_avg{2} + ka_avg{5})/2;

%% Best index
DF_dev_best_index = cell(1,5);
DF_max_dev_best_index = cell(1,5);
DF_min_dev_best_index = cell(1,5);
%Calculating max and min deviations from average within each set:
for m = 1:length(DF{1})
    for l = 1:5
        if l == 2
            DF_dev_best_index{l}(1:4,m) = [DF{5}(indexnumber_in_ka{5}(m)); ...
                DF{10}(indexnumber_in_ka{10}(m)); ...
                DF{11}(indexnumber_in_ka{11}(m)); ...
                DF{12}(indexnumber_in_ka{12}(m))];
        elseif l == 3
            DF_dev_best_index{l}(1:4,m) = [DF{6}(indexnumber_in_ka{6}(m)); ...
                DF{7}(indexnumber_in_ka{7}(m)); ...
                DF{8}(indexnumber_in_ka{8}(m)); ...
                DF{9}(indexnumber_in_ka{9}(m))];
        else
            DF_dev_best_index{l}(1:4,m) = [DF{(4*(l-1)+1)}(indexnumber_in_ka{(4*(l-1)+1)}(m)); ...
                DF{(4*(l-1)+2)}(indexnumber_in_ka{(4*(l-1)+2)}(m)); ...
                DF{(4*(l-1)+3)}(indexnumber_in_ka{(4*(l-1)+3)}(m)); ...
                DF{(4*(l-1)+4)}(indexnumber_in_ka{(4*(l-1)+4)}(m))];
        end
        DF_max_dev_best_index{l}(m) = max(DF_dev_best_index{l}(:,m),[],1);
        DF_min_dev_best_index{l}(m) = min(DF_dev_best_index{l}(:,m),[],1);
    end
end

DI_avg_best_index = cell(1,5);
DI_upper_dev_best_index = cell(1,5);
DI_lower_dev_best_index = cell(1,5);
%Turning Directivity Factor into Directivity Index:
for l = 1:5
    DI_avg_best_index{l} = 10*log10(DF_avg_best_index{l});
    DI_upper_dev_best_index{l} = 10*log10(DF_max_dev_best_index{l});
    DI_lower_dev_best_index{l} = 10*log10(DF_min_dev_best_index{l});
end

% Comparison of DI; difference

%S2 minus S1
DI_mean_diff_A_best_index = DI_avg_best_index{1} - DI_avg_best_index{2};
DI_max_diff_A_best_index = DI_upper_dev_best_index{1} - DI_lower_dev_best_index{2};

%Ventilation
DI_mean_diff_B_best_index = DI_avg_best_index{2} - DI_avg_best_index{3};
DI_max_diff_B_best_index = DI_upper_dev_best_index{2} - DI_lower_dev_best_index{3};

%Adjusted mic
DI_mean_diff_C_best_index = DI_avg_best_index{2} - DI_avg_best_index{4};
DI_max_diff_C_best_index = DI_upper_dev_best_index{2} - DI_lower_dev_best_index{4};

%Adjusted loudspeaker
DI_mean_diff_D_best_index = DI_avg_best_index{2} - DI_avg_best_index{5};
DI_max_diff_D_best_index = DI_upper_dev_best_index{2} - DI_lower_dev_best_index{5};

%% Figures

figure(3)
fig3 = semilogx(ka_avg_A,(DI_mean_diff_A-DI_mean_diff_A_best_index),'-',...
    ka_avg_A,(DI_max_diff_A-DI_max_diff_A_best_index),'-.',...
    ka_avg_B,(DI_mean_diff_B-DI_mean_diff_B_best_index),'-',...
    ka_avg_B,(DI_max_diff_B-DI_max_diff_B_best_index),'-.',...
    ka_avg_C,abs(DI_mean_diff_C-DI_mean_diff_C_best_index),'-',...
    ka_avg_C,abs(DI_max_diff_C-DI_max_diff_C_best_index),'-.',...
    ka_avg_D,abs(DI_mean_diff_D-DI_mean_diff_D_best_index),'-',...
    ka_avg_D,abs(DI_max_diff_D-DI_max_diff_D_best_index),'-.');
set(fig3,'Linewidth',1)
fig3(1).Color = [0.8500 0.3250 0.0980]; %red
fig3(2).Color = [0.8500 0.3250 0.0980]; %red
fig3(3).Color = [0 0.4470 0.7410]; %blue
fig3(4).Color = [0 0.4470 0.7410]; %blue
fig3(5).Color = [0.9290 0.6940 0.1250]; %yellow
fig3(6).Color = [0.9290 0.6940 0.1250]; %yellow
fig3(7).Color = [0.4660 0.6740 0.1880]; %green
fig3(8).Color = [0.4660 0.6740 0.1880]; %green
xlim([1.699823433909710 139.4220276990395])
xlabel({'kl, [1]',' '})
g = legend('Difference, mean A','Difference, max A',...
    'Difference, mean B','Difference, max B',...
    'Difference, mean C','Difference, max C',...
    'Difference, mean D','Difference, max D');
set(g,'Location','best');
ylabel({'Difference in directivity index, [dB]'})
set(gcf,'color','w');
grid on
title('Difference between same index and best index')


% Best index plots
figure(7)
fig6 = semilogx(...
    ka{5},DI_avg_best_index{2},'-',ka{1},DI_avg_best_index{1},'-',...
    ka{5},DI_upper_dev_best_index{2},'-.',...
    ka{5},DI_lower_dev_best_index{2},'-.',...
    ka{1},DI_upper_dev_best_index{1},'-.',...
    ka{1},DI_lower_dev_best_index{1},'-.');
set(fig6,'Linewidth',1)
fig6(1).Color = [0 0.4470 0.7410]; %blue
fig6(2).Color = [0.8500 0.3250 0.0980]; %red
fig6(3).Color = [0 0.4470 0.7410]; %blue
fig6(4).Color = [0 0.4470 0.7410]; %blue
fig6(5).Color = [0.8500 0.3250 0.0980]; %red
fig6(6).Color = [0.8500 0.3250 0.0980]; %red
% title({' ','Average Directivity Index with upper and lower',...
%     'bounds of deviation for source stationary and reinstalled',''})
xlim([1.699823433909710 139.4220276990395])
ylim([-2 8]);
xlabel({'kl, [1]'})
ylabel({'Difference in DI, [dB]'})
set(gcf,'color','w');
g = legend(...
    'LS#1',...
    'LS#2');
set(g,'Location','best');
grid on


figure(15)
subplot(4,1,3);
fig151 = semilogx(ka_avg_C,abs(DI_mean_diff_C_best_index),'-',ka_avg_C,...
    abs(DI_max_diff_C_best_index),'-.');
set(fig151,'Linewidth',1)
xlim([1.699823433909710 139.4220276990395])
ylim([-0.3 1.5])
fig151(1).Color = [0.9290 0.6940 0.1250]; %yellow
fig151(2).Color = [0.9290 0.6940 0.1250]; %yellow
% xlabel({'kl, [1]'})
% ylabel({'Difference in DI, [dB]'})
set(gcf,'color','w');
set(gca,'fontsize',10)
ylabel({'                                           Difference in directivity index, [dB]',' '})
g = legend('Mic. repositioning (mean)','Mic. repositioning (max)');
grid on

subplot(4,1,4); 
fig152 = semilogx(ka_avg_D,abs(DI_mean_diff_D_best_index),'-',...
    ka_avg_D,abs(DI_max_diff_D_best_index),'-.');
set(fig152,'Linewidth',1)
xlim([1.699823433909710 139.4220276990395]) %1.980014109828893
ylim([-0.3 1.5])
fig152(1).Color = [0.4660 0.6740 0.1880]; %green
fig152(2).Color = [0.4660 0.6740 0.1880]; %green
% xlabel({'kl, [1]'})
% ylabel({'Difference in DI, [dB]'})
set(gcf,'color','w');
set(gca,'fontsize',10)
xlabel({'kl, [1]'})
g = legend('LS repositioning (mean)', 'LS repositioning (max)');
%set(g,'Location','best');
grid on

subplot(4,1,1);
fig153 = semilogx(ka_avg_A,DI_mean_diff_A_best_index,'-',ka_avg_A,...
    DI_max_diff_A_best_index,'-.');
set(fig153,'Linewidth',1)
xlim([1.699823433909710 139.4220276990395])
ylim([-0.3 1.5])
fig153(1).Color = [0 0.4470 0.7410]; %blue
fig153(2).Color = [0 0.4470 0.7410];
%xlabel({'kl, [1]'})
set(gcf,'color','w');
set(gca,'fontsize',10)
g = legend('LS#2 vs. LS#1 (mean)','LS#2 vs. LS#1 (max)');
grid on

subplot(4,1,2); 
fig154 = semilogx(ka_avg_B,DI_mean_diff_B_best_index,'-',ka_avg_B,...
    DI_max_diff_B_best_index,'-.');
set(fig154,'Linewidth',1)
xlim([1.699823433909710 139.4220276990395])
ylim([-0.3 1.5])
fig154(1).Color = [0.8500 0.3250 0.0980]; %red
fig154(2).Color = [0.8500 0.3250 0.0980];
% ylabel({'Difference in DI, [dB]'})
set(gcf,'color','w');
set(gca,'fontsize',10)
g = legend('Vent. on vs. off (mean)','Vent. on vs. off (max)');
grid on



% Plots before editing for best index:
% figure(7)
% fig6 = semilogx(...
%     ka{5},DI_avg{2},'-',ka{1},DI_avg{1},'-',...
%     ka{5},DI_upper_dev{2},'-.',...
%     ka{5},DI_lower_dev{2},'-.',...
%     ka{1},DI_upper_dev{1},'-.',...
%     ka{1},DI_lower_dev{1},'-.');
% set(fig6,'Linewidth',1)
% fig6(1).Color = [0 0.4470 0.7410]; %blue
% fig6(2).Color = [0.8500 0.3250 0.0980]; %red
% fig6(3).Color = [0 0.4470 0.7410]; %blue
% fig6(4).Color = [0 0.4470 0.7410]; %blue
% fig6(5).Color = [0.8500 0.3250 0.0980]; %red
% fig6(6).Color = [0.8500 0.3250 0.0980]; %red
% % title({' ','Average Directivity Index with upper and lower',...
% %     'bounds of deviation for source stationary and reinstalled',''})
% xlim([1.699823433909710 139.4220276990395])
% ylim([-2 8]);
% xlabel({'kl, [1]'})
% ylabel({'Difference in DI, [dB]'})
% set(gcf,'color','w');
% g = legend(...
%     'LS#1',...
%     'LS#2');
% set(g,'Location','best');
% grid on

% 
% figure(15)
% subplot(4,1,3);
% fig151 = semilogx(ka_avg_C,abs(DI_mean_diff_C),'-',ka_avg_C,...
%     abs(DI_max_diff_C),'-.');
% set(fig151,'Linewidth',1)
% xlim([1.699823433909710 139.4220276990395])
% ylim([-0.3 1.5])
% fig151(1).Color = [0.9290 0.6940 0.1250]; %yellow
% fig151(2).Color = [0.9290 0.6940 0.1250]; %yellow
% % xlabel({'kl, [1]'})
% % ylabel({'Difference in DI, [dB]'})
% set(gcf,'color','w');
% set(gca,'fontsize',10)
% ylabel({'                                           Difference in directivity index, [dB]',' '})
% g = legend('Mic. repositioning (mean)','Mic. repositioning (max)');
% grid on
% 
% subplot(4,1,4); 
% fig152 = semilogx(ka_avg_D,abs(DI_mean_diff_D),'-',...
%     ka_avg_D,abs(DI_max_diff_D),'-.');
% set(fig152,'Linewidth',1)
% xlim([1.699823433909710 139.4220276990395]) %1.980014109828893
% ylim([-0.3 1.5])
% fig152(1).Color = [0.4660 0.6740 0.1880]; %green
% fig152(2).Color = [0.4660 0.6740 0.1880]; %green
% % xlabel({'kl, [1]'})
% % ylabel({'Difference in DI, [dB]'})
% set(gcf,'color','w');
% set(gca,'fontsize',10)
% xlabel({'kl, [1]'})
% g = legend('LS repositioning (mean)', 'LS repositioning (max)');
% %set(g,'Location','best');
% grid on
% 
% subplot(4,1,1);
% fig153 = semilogx(ka_avg_A,DI_mean_diff_A,'-',ka_avg_A,DI_max_diff_A,...
%     '-.');
% set(fig153,'Linewidth',1)
% xlim([1.699823433909710 139.4220276990395])
% ylim([-0.3 1.5])
% fig153(1).Color = [0 0.4470 0.7410]; %blue
% fig153(2).Color = [0 0.4470 0.7410];
% %xlabel({'kl, [1]'})
% set(gcf,'color','w');
% set(gca,'fontsize',10)
% g = legend('LS#1 vs. LS#2 (mean)','LS#1 vs. LS#2 (max)');
% grid on
% 
% subplot(4,1,2); 
% fig154 = semilogx(ka_avg_B,DI_mean_diff_B,'-',ka_avg_B,DI_max_diff_B,'-.');
% set(fig154,'Linewidth',1)
% xlim([1.699823433909710 139.4220276990395])
% ylim([-0.3 1.5])
% fig154(1).Color = [0.8500 0.3250 0.0980]; %red
% fig154(2).Color = [0.8500 0.3250 0.0980];
% % ylabel({'Difference in DI, [dB]'})
% set(gcf,'color','w');
% set(gca,'fontsize',10)
% g = legend('Vent. on vs. off (mean)','Vent. on vs. off (max)');
% grid on
% 



% figure(13)
% subplot(2,1,1);
% fig131 = semilogx(ka_avg{2},abs(DI_mean_diff_C),'-',ka_avg{2},...
%     abs(DI_max_diff_C),'-.');
% %title({' ','Difference in Directivity Indices, absolute values',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig131,'Linewidth',1)
% %xlim([1.699823433909710 180])
% xlim([1.699823433909710 139.4220276990395])
% ylim([-0.3 1.5])
% fig131(1).Color = [0.9290 0.6940 0.1250]; %yellow
% fig131(2).Color = [0.9290 0.6940 0.1250]; %yellow
% xlabel({'kl, [1]',' '})
% ylabel({' ','Difference in DI, [dB]',' '})
% set(gcf,'color','w');
% g = legend('Mic. repositioning (mean)','Mic. repositioning (max)');
% %set(g,'Location','best');
% grid on
% 
% subplot(2,1,2); 
% fig132 = semilogx(ka_avg{2},abs(DI_mean_diff_D),'-',...
%     ka_avg{2},abs(DI_max_diff_D),'-.');
% %title({' ','Difference in Directivity Indices, absolute values',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig132,'Linewidth',1)
% %xlim([1.699823433909710 180])
% xlim([1.699823433909710 139.4220276990395]) %1.980014109828893
% ylim([-0.3 1.5])
% fig132(1).Color = [0.4660 0.6740 0.1880]; %green
% fig132(2).Color = [0.4660 0.6740 0.1880]; %green
% xlabel({'kl, [1]',' '})
% ylabel({' ','Difference in DI, [dB]',' '})
% set(gcf,'color','w');
% g = legend('LS repositioning (mean)', 'LS repositioning (max)');
% %set(g,'Location','best');
% grid on

% figure(14)
% subplot(2,1,1);
% fig141 = semilogx(ka_avg{2},DI_mean_diff_A,'-',ka_avg{2},DI_max_diff_A,...
%     '-.');
% %title({' ','Difference in Directivity Indices, absolute values',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig141,'Linewidth',1)
% %xlim([1.699823433909710 180])
% xlim([1.699823433909710 139.4220276990395])
% ylim([-0.3 1.5])
% fig141(1).Color = [0 0.4470 0.7410]; %blue
% fig141(2).Color = [0 0.4470 0.7410];
% xlabel({'kl, [1]'})
% ylabel({'Difference in DI, [dB]'})
% set(gcf,'color','w');
% g = legend('LS#1 vs. LS#2 (mean)','LS#1 vs. LS#2 (max)');
% grid on
% 
% subplot(2,1,2); 
% fig142 = semilogx(ka_avg{2},DI_mean_diff_B,'-',ka_avg{2},DI_max_diff_B,'-.');
% %title({' ','Difference in Directivity Indices, absolute values',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig142,'Linewidth',1)
% %xlim([1.699823433909710 180])
% xlim([1.699823433909710 139.4220276990395])
% ylim([-0.3 1.5])
% fig142(1).Color = [0.8500 0.3250 0.0980]; %red
% fig142(2).Color = [0.8500 0.3250 0.0980];
% xlabel({'kl, [1]'})
% ylabel({'Difference in DI, [dB]'})
% set(gcf,'color','w');
% g = legend('Vent. on vs. off (mean)','Vent. on vs. off (max)');
% %set(g,'Location','best');
% grid on


%% Plots

% figure(8)
% fig8 = semilogx(ka_avg{2},DI_mean_diff_A,'-',ka_avg{2},DI_max_diff_A,...
%     '-.',ka_avg{2},DI_mean_diff_B,'-',ka_avg{2},DI_max_diff_B,'-.',...
%     ka_avg{2},DI_mean_diff_C,'-',ka_avg{2},DI_max_diff_C,'-.',...
%     ka_avg{2},DI_mean_diff_D,'-',ka_avg{2},DI_max_diff_D,'-.');
% title({' ','Difference in Directivity Indices',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig8,'Linewidth',1)
% xlim([1.3 100])
% fig8(1).Color = [0 0.4470 0.7410]; %blue
% fig8(2).Color = [0 0.4470 0.7410];
% fig8(3).Color = [0.8500 0.3250 0.0980]; %red
% fig8(4).Color = [0.8500 0.3250 0.0980];
% fig8(5).Color = [0.9290 0.6940 0.1250]; %yellow
% fig8(6).Color = [0.9290 0.6940 0.1250]; %yellow
% fig8(7).Color = [0.4660 0.6740 0.1880]; %green
% fig8(8).Color = [0.4660 0.6740 0.1880]; %green
% xlabel({'kl, [1]',' '})
% ylabel({' ','Difference, [dB]',' '})
% set(gcf,'color','w');
% g = legend('LS#1 vs. LS#2 (mean)','LS#1 vs. LS#2 (max)',...
%     'Vent. on vs. off (mean)','Vent. on vs. off (max)',...
%     'Mic. adjust (mean)','Mic. adjusted (max)',...
%     'Loudspeaker adjusted (mean)', 'Loudspeaker adjusted (max)');
% set(g,'Location','best');
% grid on
% 
% figure(9)
% fig9 = semilogx(ka_avg{2},DI_mean_diff_A,'-',ka_avg{2},DI_max_diff_A,...
%     '-.',ka_avg{2},DI_mean_diff_B,'-',ka_avg{2},DI_max_diff_B,'-.');
% title({' ','Difference in Directivity Indices',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig9,'Linewidth',1)
% xlim([1.3 100])
% fig9(1).Color = [0 0.4470 0.7410]; %blue
% fig9(2).Color = [0 0.4470 0.7410];
% fig9(3).Color = [0.8500 0.3250 0.0980]; %red
% fig9(4).Color = [0.8500 0.3250 0.0980];
% xlabel({'kl, [1]',' '})
% ylabel({' ','Difference, [dB]',' '})
% set(gcf,'color','w');
% g = legend('LS#1 vs. LS#2 (mean)','LS#1 vs. LS#2 (max)',...
%     'Vent. on vs. off (mean)','Vent. on vs. off (max)');
% set(g,'Location','best');
% grid on
% 
% figure(10)
% fig10 = semilogx(ka_avg{2},abs(DI_mean_diff_C),'-',ka_avg{2},...
%     abs(DI_max_diff_C),'-.',ka_avg{2},abs(DI_mean_diff_D),'-',...
%     ka_avg{2},abs(DI_max_diff_D),'-.');
% title({' ','Difference in Directivity Indices, absolute values',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig10,'Linewidth',1)
% xlim([1.3 100])
% ylim([-0.2 1.2])
% fig10(1).Color = [0.9290 0.6940 0.1250]; %yellow
% fig10(2).Color = [0.9290 0.6940 0.1250]; %yellow
% fig10(3).Color = [0.4660 0.6740 0.1880]; %green
% fig10(4).Color = [0.4660 0.6740 0.1880]; %green
% xlabel({'kl, [1]',' '})
% ylabel({' ','Difference, [dB]',' '})
% set(gcf,'color','w');
% g = legend('Mic. adjust (mean)','Mic. adjusted (max)',...
%     'Loudspeaker adjusted (mean)', 'Loudspeaker adjusted (max)');
% set(g,'Location','best');
% grid on
% 
% figure(12)
% fig12 = semilogx(ka_avg{2},abs(DI_mean_diff_C),'-',ka_avg{2},...
%     abs(DI_max_diff_C),'-.',ka_avg{2},abs(DI_mean_diff_D),'-',...
%     ka_avg{2},abs(DI_max_diff_D),'-.');
% title({' ','Difference in Directivity Indices, absolute values',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig12,'Linewidth',1)
% % xlim([1.3 100])
% % ylim([-0.2 1.2])
% fig12(1).Color = [0.9290 0.6940 0.1250]; %yellow
% fig12(2).Color = [0.9290 0.6940 0.1250]; %yellow
% fig12(3).Color = [0.4660 0.6740 0.1880]; %green
% fig12(4).Color = [0.4660 0.6740 0.1880]; %green
% xlabel({'kl, [1]',' '})
% ylabel({' ','Difference, [dB]',' '})
% set(gcf,'color','w');
% g = legend('Mic. adjust (mean)','Mic. adjusted (max)',...
%     'Loudspeaker adjusted (mean)', 'Loudspeaker adjusted (max)');
% set(g,'Location','best');
% grid on
% 
% figure(11)
% fig11 = semilogx(ka_avg{2},(DI_mean_diff_C),'-',ka_avg{2},...
%     (DI_max_diff_C),'-.',ka_avg{2},(DI_mean_diff_D),'-',...
%     ka_avg{2},(DI_max_diff_D),'-.');
% title({' ','Difference in Directivity Indices',' '})
% %    'between loudspeaker #1 & #2',''})
% set(fig11,'Linewidth',1)
% xlim([1.3 100])
% fig11(1).Color = [0.9290 0.6940 0.1250]; %yellow
% fig11(2).Color = [0.9290 0.6940 0.1250]; %yellow
% fig11(3).Color = [0.4660 0.6740 0.1880]; %green
% fig11(4).Color = [0.4660 0.6740 0.1880]; %green
% xlabel({'kl, [1]',' '})
% ylabel({' ','Difference, [dB]',' '})
% set(gcf,'color','w');
% g = legend('Mic. adjust (mean)','Mic. adjusted (max)',...
%     'Loudspeaker adjusted (mean)', 'Loudspeaker adjusted (max)');
% set(g,'Location','best');
% grid on


% %Plotting Directivity Index for measurement setup A and B with upper and
% %lower bounds of deviation:
% figure(1)
% semilogx(ka_avg{2},DI_avg{2},'-',ka_avg{2},DI_upper_dev{2},...
%     '-.',ka_avg{2},DI_lower_dev{2},'-.',...
%     ka_avg{1},DI_avg{1},'-',ka_avg{1},DI_upper_dev{1},'-.',ka_avg{1},...
%     DI_lower_dev{1},'-.');
% title({' ','Average Directivity Index with upper and lower',...
%     'bounds of deviation for loudspeaker #1 & #2',''})
% xlabel({'ka, [1]',' '})
% ylabel({' ','Magnitude, dB',' '})
% set(gcf,'color','w');
% g = legend('LS #1','Upper dev.','Lower dev.','LS #2',...
%     'Upper dev.','Lower dev.');
% set(g,'Location','best');
% grid on

% % % % % % figure(2)
% % % % % % semilogx(ka{1},DI{1},'-o',ka{2},DI{2},'-o',ka{3},DI{3},'-o',ka{4},DI{4},'-o',...
% % % % % %     ka_avg{1},DI_avg{1},'-o')%,ka_avg{1},(DI_upper_dev{1}),'-.',...
% % % % % %     %ka_avg{1},(DI_lower_dev{1}),'-.')
% % % % % % xlabel({'ka, [1]',' '})
% % % % % % ylabel({' ','Directivity Index, [dB]',' '})
% % % % % % legend('Series 1','Series 2','Series 3', 'Series 4','Average')%,'Max dev.',...
% % % % % %     %'Min dev.');
% % % % % % set(gcf,'color','w');
% % % % % % title({' ','Average Directivity Factor with upper and lower',...
% % % % % %     'bounds of deviation for loudspeaker #2',''})
% % % % % % grid on
% % % % % 
% % % % % % figure(7)
% % % % % % semilogx(ka_avg{1},DI{1},'-o',ka_avg{1},DI{2},'-o',ka_avg{1},DI{3},'-o',ka_avg{1},DI{4},'-o',...
% % % % % %     ka_avg{1},DI_avg{1},'-o')%,ka_avg{1},(DI_upper_dev{1}),'-.',...
% % % % % %     %ka_avg{1},(DI_lower_dev{1}),'-.')
% % % % % % xlabel({'ka, [1]',' '})
% % % % % % ylabel({' ','Directivity Index, [dB]',' '})
% % % % % % legend('Series 1','Series 2','Series 3', 'Series 4','Average')%,'Max dev.',...
% % % % % %     %'Min dev.');
% % % % % % set(gcf,'color','w');
% % % % % % title({' ','Average Directivity Factor with upper and lower',...
% % % % % %     'bounds of deviation for loudspeaker #2',''})
% % % % % % grid on

% figure(2)
% loglog(ka{1},DF{1},'-o',ka{2},DF{2},'-o',ka{3},DF{3},'-o',ka{4},DF{4},'-o',...
%     ka_avg{1},DF_avg{1},'-o')%,ka_avg{1},(DF_max_dev{1}),'-.',...
%     %ka_avg{1},(DF_min_dev{1}),'-.')
% xlabel({'ka, [1]',' '})
% ylabel({' ','Directivity Factor, [1]',' '})
% legend('Series 1','Series 2','Series 3', 'Series 4','Average')%,'Max dev.',...
%     %'Min dev.');
% set(gcf,'color','w');
% title({' ','Directivity Factor for each measurement series performed on',...
%     'box for loudspeaker #2, including average Directivity Factor',''})
% grid on

% figure(7)
% loglog(ka_avg{1},DF{1},'-o',ka_avg{1},DF{2},'-o',ka_avg{1},DF{3},'-o',ka_avg{1},DF{4},'-o',...
%     ka_avg{1},DF_avg{1},'-o')%,ka_avg{1},(DF_max_dev{1}),'-.',...
%     %ka_avg{1},(DF_min_dev{1}),'-.')
% xlabel({'ka, [1]',' '})
% ylabel({' ','Directivity Factor, [1]',' '})
% legend('Series 1','Series 2','Series 3', 'Series 4','Average')%,'Max dev.',...
%     %'Min dev.');
% set(gcf,'color','w');
% title({' ','Directivity Factor for each measurement series performed on',...
%     'box for loudspeaker #2, including average Directivity Factor',...
%     'by use of average ka numbers for points in series',''})
% grid on


% figure(3)
% fig3 = semilogx(...
%     ka{5},DI_avg{2},'-',ka{1},DI_avg{1},'-',ka{6},DI_avg{3},'-',...
%     ka{13},DI_avg{4},'-',ka{17},DI_avg{5},'-',...
%     ka{5},DI_upper_dev{2},'-.',...
%     ka{5},DI_lower_dev{2},'-.',...
%     ka{1},DI_upper_dev{1},'-.',...
%     ka{1},DI_lower_dev{1},'-.',...
%     ka{6},DI_upper_dev{3},'-.',...
%     ka{6},DI_lower_dev{3},'-.',...
%     ka{13},DI_upper_dev{4},'-.',...
%     ka{13},DI_lower_dev{4},'-.',...
%     ka{17},DI_upper_dev{5},'-.',...
%     ka{17},DI_lower_dev{5},'-.');
% set(fig3,'Linewidth',1)
% fig3(1).Color = [0 0.4470 0.7410]; %blue
% fig3(6).Color = [0 0.4470 0.7410]; %blue
% fig3(7).Color = [0 0.4470 0.7410]; %blue
% fig3(2).Color = [0.8500 0.3250 0.0980]; %red
% fig3(8).Color = [0.8500 0.3250 0.0980]; %red
% fig3(9).Color = [0.8500 0.3250 0.0980]; %red
% fig3(3).Color = [0.9290 0.6940 0.1250]; %yellow
% fig3(10).Color = [0.9290 0.6940 0.1250]; %yellow
% fig3(11).Color = [0.9290 0.6940 0.1250]; %yellow
% fig3(4).Color = [0.4660 0.6740 0.1880]; %green
% fig3(12).Color = [0.4660 0.6740 0.1880]; %green
% fig3(13).Color = [0.4660 0.6740 0.1880]; %green
% fig3(5).Color = [0.4940 0.1840 0.5560]; %purple
% fig3(14).Color = [0.4940 0.1840 0.5560]; %purple
% fig3(15).Color = [0.4940 0.1840 0.5560]; %purple
% % title({' ','Average Directivity Index with upper and lower',...
% %     'bounds of deviation for setups 1-5 for box',''})
% xlabel({'ka, [1]'})
% ylabel({'Directivity Index, [dB]'})
% set(gcf,'color','w');
% g = legend('LS#1','LS#2','LS#1 vent. on','LS#1 R reinstalled',...
%     'LS#1 S reinstalled');
% set(g,'Location','best');
% grid on

% figure(4)
% semilogx(...
%     ka{5},DI_avg{2},'-',ka{5},DI_upper_dev{2},'-.',...
%     ka{5},DI_lower_dev{2},'-.',...
%     ka{6},DI_avg{3},'-',ka{6},DI_upper_dev{3},'-.',...
%     ka{6},DI_lower_dev{3},'-.');
% title({' ','Average Directivity Index with upper and lower',...
%     'bounds of deviation for ventillation on and off',''})
% xlabel({'ka, [1]',' '})
% ylabel({' ','Directivity Index, dB',' '})
% set(gcf,'color','w');
% g = legend(...
%     'LS #1 vent. off','Upper dev.','Lower dev.',...
%     'LS #1 vent. on','Upper dev.','Lower dev.');
% set(g,'Location','best');
% grid on

% figure(5)
% semilogx(...
%     ka{5},DI_avg{2},'-',ka{5},DI_upper_dev{2},'-.',...
%     ka{5},DI_lower_dev{2},'-.',...
%     ka{13},DI_avg{4},'-',ka{13},DI_upper_dev{4},'-.',...
%     ka{13},DI_lower_dev{4},'-.');
% title({' ','Average Directivity Index with upper and lower',...
%     'bounds of deviation for box stationary and reinstalled',''})
% xlabel({'ka, [1]',' '})
% ylabel({' ','Directivity Index, dB',' '})
% set(gcf,'color','w');
% g = legend(...
%     'LS #1, R stationary','Upper dev.','Lower dev.',...
%     'LS #1, R reinstalled','Upper dev.','Lower dev.');
% set(g,'Location','best');
% grid on

% figure(6)
% fig6 = semilogx(...
%     ka{5},DI_avg{2},'-',ka{17},DI_avg{5},'-',...
%     ka{5},DI_upper_dev{2},'-.',...
%     ka{5},DI_lower_dev{2},'-.',...
%     ka{17},DI_upper_dev{5},'-.',...
%     ka{17},DI_lower_dev{5},'-.');
% set(fig6,'Linewidth',1)
% fig6(1).Color = [0 0.4470 0.7410]; %blue
% fig6(2).Color = [0.8500 0.3250 0.0980]; %red
% fig6(3).Color = [0 0.4470 0.7410]; %blue
% fig6(4).Color = [0 0.4470 0.7410]; %blue
% fig6(5).Color = [0.8500 0.3250 0.0980]; %red
% fig6(6).Color = [0.8500 0.3250 0.0980]; %red
% % title({' ','Average Directivity Index with upper and lower',...
% %     'bounds of deviation for source stationary and reinstalled',''})
% xlim([1.980014109828893 250]); %Limiting from ka value cutoff from SNR
% ylim([-9 13]);
% xlabel({'ka, [1]'})
% ylabel({'Directivity index, [dB]'})
% set(gcf,'color','w');
% g = legend(...
%     'LS#1, S stationary',...
%     'LS#1, S reinstalled');
% set(g,'Location','best');
% grid on

