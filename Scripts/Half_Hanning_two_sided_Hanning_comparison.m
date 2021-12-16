%Script to test whether adding a fade-in to the half-Hanning window, making
%it a two-sided window, will change the results for the Directivity Index.
%Tested only on data for the cylinder, as the box data would require a
%shift in the window placement, making the results less comparable.

%Written 30.11.21 by Isabel Berg

%FFTset_v4_two_sided_Hanning (taper in)
%FFTset_v4

close all
clear all

%% Comparison for cylinder measurements (1):
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

%S1 and S2
DI_mean_diff_A = DI_avg{2} - DI_avg{1};
DI_max_diff_A = DI_upper_dev{2} - DI_lower_dev{1};

%% Comparison for cylinder measurements (2): (with taper in)

%Perfoming Fourier Transformation on Impulse Responses:
for k = 1:8
    if k>=4
        lspnum = 2;
    else
        lspnum = 1;
    end
    [TF_ti{k},DF_nom_ti{k},DF_denom_ti{k},DI_ti{k},ka_ti{k},fvec_ka_ti{k}]= ...
        FFTset_v4_two_sided_Hanning(dataset{k},lspnum,cdata_avg(k,:),radius);
end

%Calculating Directivity Factor:
DF_ti = cell(1,numfiles);
for k = 1:numfiles
    DF_ti{k} = DF_nom_ti{k}./DF_denom_ti{k};
end

DF_avg_ti = cell(1,2);
ka_avg_ti = cell(1,2);
%Averaging over four series for each setup:
for l = 1:2
    if l == 1
        DF_avg_ti{l} = (DF_ti{1} + DF_ti{2} + DF_ti{3} + DF_ti{4})/4;
        ka_avg_ti{l} = (ka_ti{1} + ka_ti{2} + ka_ti{3} + ka_ti{4})/4;
    else
        DF_avg_ti{l} = (DF_ti{5} + DF_ti{6} + DF_ti{7} + DF_ti{8})/4;
        ka_avg_ti{l} = (ka_ti{5} + ka_ti{6} + ka_ti{7} + ka_ti{8})/4;
    end
end

DF_dev_ti = cell(1,2);
DF_max_dev_ti = cell(1,2);
DF_min_dev_ti = cell(1,2);
%Calculating max and min deviations from average within each set:
for l = 1:2
    if l == 1
        DF_dev_ti{l} = [(DF_ti{1}); DF_ti{2}; DF_ti{3}; DF_ti{4}];
    else
        DF_dev_ti{l} = [(DF_ti{5}); DF_ti{6}; DF_ti{7}; DF_ti{8}];
    end
    DF_max_dev_ti{l} = max(DF_dev_ti{l},[],1);
    DF_min_dev_ti{l} = min(DF_dev_ti{l},[],1);
end

DI_avg_ti = cell(1,2);
DI_upper_dev_ti = cell(1,2);
DI_lower_dev_ti = cell(1,2);
%Turning Directivity Factor into Directivity Index:
for l = 1:2
    DI_avg_ti{l} = 10*log10(DF_avg{l});
    DI_upper_dev_ti{l} = 10*log10(DF_max_dev{l});
    DI_lower_dev_ti{l} = 10*log10(DF_min_dev{l});
end

%S1 and S2
DI_mean_diff_A_ti = DI_avg_ti{2} - DI_avg_ti{1};
DI_max_diff_A_ti = DI_upper_dev_ti{2} - DI_lower_dev_ti{1};


diffka = ka_avg{1} - ka_avg_ti{1};
diffDImeandiff = DI_mean_diff_A - DI_mean_diff_A_ti;


%% Figures

% figure(14)
% fig141 = semilogx(ka_avg{1},DI_mean_diff_A,'-o',ka_avg{1},DI_max_diff_A,...
%     '-.',ka_avg_ti{1},DI_mean_diff_A_ti,'-o',ka_avg_ti{1},DI_max_diff_A_ti,...
%     '-.');
% set(fig141,'Linewidth',1)
% xlim([0.435441101586095 43.853587469566044])
% ylim([-0.3 1.5])
% fig141(1).Color = [0.9290 0.6940 0.1250]; %yellow
% fig141(2).Color = [0.9290 0.6940 0.1250];
% fig141(3).Color = [0.4660 0.6740 0.1880]; %green
% fig141(4).Color = [0.4660 0.6740 0.1880];
% xlabel({'ka, [1]',' '})
% ylabel({'Difference in directivity index, [dB]'})
% set(gcf,'color','w');
% g = legend('LS#1 vs. LS#2 (mean)','LS#1 vs. LS#2 (max)',...
%     'LS#1 vs. LS#2 (mean), taper in','LS#1 vs. LS#2 (max), taper in');
% grid on
% 
% 
% figure(7)
% fig6 = semilogx(ka{1},DI_avg{1},'-o',...
%     ka{5},DI_avg{2},'-o',...
%     ka_ti{1},DI_avg_ti{1},'-o',...
%     ka_ti{5},DI_avg_ti{2},'-o',...
%     ka{1},DI_upper_dev{1},'-.',...
%     ka{1},DI_lower_dev{1},'-.',...
%     ka{5},DI_upper_dev{2},'-.',...
%     ka{5},DI_lower_dev{2},'-.',...
%     ka_ti{1},DI_upper_dev_ti{1},'-.',...
%     ka_ti{1},DI_lower_dev_ti{1},'-.',...
%     ka_ti{5},DI_upper_dev_ti{2},'-.',...
%     ka_ti{5},DI_lower_dev_ti{2},'-.');
% set(fig6,'Linewidth',1)
% fig6(1).Color = [0 0.4470 0.7410]; %blue
% fig6(2).Color = [0.8500 0.3250 0.0980]; %red
% fig6(3).Color = [0.9290 0.6940 0.1250]; %yellow
% fig6(4).Color = [0.4660 0.6740 0.1880]; %green
% fig6(5).Color = [0 0.4470 0.7410]; %blue
% fig6(6).Color = [0 0.4470 0.7410]; %blue
% fig6(7).Color = [0.8500 0.3250 0.0980]; %red
% fig6(8).Color = [0.8500 0.3250 0.0980]; %red
% fig6(9).Color = [0.9290 0.6940 0.1250]; %yellow
% fig6(10).Color = [0.9290 0.6940 0.1250]; %yellow
% fig6(11).Color = [0.4660 0.6740 0.1880]; %green
% fig6(12).Color = [0.4660 0.6740 0.1880]; %green
% xlim([0.435441101586095 43.853587469566044])
% ylim([1 8]);
% xlabel({'ka, [1]'})
% ylabel({'Directivity index, [dB]'})
% set(gcf,'color','w');
% g = legend(...
%     'LS#1',...
%     'LS#2',...
%     'LS#1 (taper in)',...
%     'LS#2 (taper in)');
% set(g,'Location','best');
% grid on
% % 
% 
% %% Lower cutoff from SNR = 30 dB:
% %LS1: 0.435441101586095 ca. 237 Hz
% %LS2: 0.204445629228746 ca. 111 Hz
% %LS2: 43.853587469566044 ca. 23879.88281250000 Hz upper limit
% 
figure(2)
plot(DI_mean_diff_A,diffka,'-o')
grid on

figure(3)
semilogx(ka_avg{1},diffDImeandiff,'-o')
grid on

clear all
%% For box:
%% Comparison for box measurements (1):
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
    [TF{k},DF_nom{k},DF_denom{k},DI{k},ka{k},fvec{k}]= ...
        FFTset_v4(dataset{k},lspnum,cdata_avg(k,:),radius);
end


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

%S1 and S2
DI_mean_diff_A = DI_avg{1} - DI_avg{2};
DI_max_diff_A = DI_upper_dev{1} - DI_lower_dev{2};

%Ventilation
DI_mean_diff_B = DI_avg{2} - DI_avg{3};
DI_max_diff_B = DI_upper_dev{2} - DI_lower_dev{3};

%Adjusted mic
DI_mean_diff_C = DI_avg{2} - DI_avg{4};
DI_max_diff_C = DI_upper_dev{2} - DI_lower_dev{4};

%Adjusted loudspeaker
DI_mean_diff_D = DI_avg{2} - DI_avg{5};
DI_max_diff_D = DI_upper_dev{2} - DI_lower_dev{5};

%% Comparison for box measurements (2):

%Perfoming Fourier Transformation on Impulse Responses:
for k = 1:20
    if k<=4
        lspnum = 2;
    else
        lspnum = 1;
    end
    [TF_ti{k},DF_nom_ti{k},DF_denom_ti{k},DI_ti{k},ka_ti{k},fvec_ti{k}]= ...
        FFTset_v4(dataset{k},lspnum,cdata_avg(k,:),radius);
end


%Calculating Directivity Factor:
DF_ti = cell(1,20);
for k = 1:20
    DF_ti{k} = DF_nom_ti{k}./DF_denom_ti{k};
end

DF_avg_ti = cell(1,5);
ka_avg_ti = cell(1,5);
%Averaging over four series for each setup:
for l = 1:5
    if l == 2
        DF_avg_ti{l} = (DF_ti{5} + DF_ti{10} + DF_ti{11} + DF_ti{12})/4;
        ka_avg_ti{l} = (ka_ti{5} + ka_ti{10} + ka_ti{11} + ka_ti{12})/4;
    elseif l == 3
        DF_avg_ti{l} = (DF_ti{6} + DF_ti{7} + DF_ti{8} + DF_ti{9})/4;
        ka_avg_ti{l} = (ka_ti{6} + ka_ti{7} + ka_ti{8} + ka_ti{9})/4;
    else
        DF_avg_ti{l} = (DF_ti{4*(l-1)+1} + DF_ti{4*(l-1)+2} + DF_ti{4*(l-1)+3} + ...
            DF_ti{4*(l-1)+4})/4;
        ka_avg_ti{l} = (ka{4*(l-1)+1} + ka{4*(l-1)+2} + ka{4*(l-1)+3} + ...
            ka{4*(l-1)+4})/4;
    end
end

DF_dev_ti = cell(1,5);
DF_max_dev_ti = cell(1,5);
DF_min_dev_ti = cell(1,5);
%Calculating max and min deviations from average within each set:
for l = 1:5
    if l == 2
        DF_dev_ti{l} = [(DF{5}); DF{10}; DF{11}; DF{12}];
    elseif l == 3
        DF_dev_ti{l} = [(DF{6}); DF{7}; DF{8}; DF{9}];
    else
        DF_dev_ti{l} = [(DF{4*(l-1)+1}); DF{4*(l-1)+2}; DF{4*(l-1)+3}; ...
            (DF{4*(l-1)+4})]; 
    end
    DF_max_dev_ti{l} = max(DF_dev_ti{l},[],1);
    DF_min_dev_ti{l} = min(DF_dev_ti{l},[],1);
end

DI_avg_ti = cell(1,5);
DI_upper_dev_ti = cell(1,5);
DI_lower_dev_ti = cell(1,5);
%Turning Directivity Factor into Directivity Index:
for l = 1:5
    DI_avg_ti{l} = 10*log10(DF_avg{l});
    DI_upper_dev_ti{l} = 10*log10(DF_max_dev{l});
    DI_lower_dev_ti{l} = 10*log10(DF_min_dev{l});
end

%S1 and S2
DI_mean_diff_A_ti = DI_avg_ti{1} - DI_avg_ti{2};
DI_max_diff_A_ti = DI_upper_dev_ti{1} - DI_lower_dev_ti{2};

%Ventilation
DI_mean_diff_B_ti = DI_avg_ti{2} - DI_avg_ti{3};
DI_max_diff_B_ti = DI_upper_dev_ti{2} - DI_lower_dev_ti{3};

%Adjusted mic
DI_mean_diff_C_ti = DI_avg_ti{2} - DI_avg_ti{4};
DI_max_diff_C_ti = DI_upper_dev_ti{2} - DI_lower_dev_ti{4};

%Adjusted loudspeaker
DI_mean_diff_D_ti = DI_avg_ti{2} - DI_avg_ti{5};
DI_max_diff_D_ti = DI_upper_dev_ti{2} - DI_lower_dev_ti{5};


diffka = ka_avg{1} - ka_avg_ti{1};
diffDImeandiff = DI_mean_diff_A - DI_mean_diff_A_ti;


%% Figures

figure(5)
plot(DI_mean_diff_A,diffka,'-o')
grid on

figure(6)
semilogx(ka_avg{1},diffDImeandiff,'-o')
grid on