% Script to calculate the simulated directivity index for different values
% of ngauss as input to the EDtoolbox. Gives one plot for 10 kHz and one
% for 20 kHz, showing when the calculated DI falls within 1/1000 of the
% value for the reference of ngauss = 68.

% Written 11.11.21 by Isabel Berg

recalc = 2;
%recalc = 1;
%recalc = 0;

frequency = [10000 20000];
ngaussvec = [5 10 15 20 25 30 35 40 42 45 46 47 48 49 50 55 60 65 68 70];
%ngaussvec = [5 10 15];    %Test ngauss values
 %% Calculations
if recalc == 1
    P_matrix = cell(length(frequency),1); %Using EDtoolbox to calculate transfer function
    for kk = 1:length(ngaussvec)
        [tftot] = (boxsimDIngauss(frequency,ngaussvec(kk)));
        for jj = 1:length(frequency) % Looping for each frequency value
            P_matrix{jj}(kk,:) = tftot(jj,:);
        end
    end

    weightvec = ones(length(tftot),1); %Weighting since only sim. 1:pi
    weigthvec(1) = 0.5;
    weightvec(end) = 0.5;
    
    DF_num_matrix = cell(length(frequency),1);
    DF_denom_matrix = cell(length(frequency),1);
    DF_matrix = cell(length(frequency),1);
    DI_matrix = cell(length(frequency),1);
    for ll = 1:length(ngaussvec)
        for jj = 1:length(frequency)
            DF_num_matrix{jj}(ll,:) = abs((P_matrix{jj}(ll,1)).^2);
            y = sum((abs((P_matrix{jj}(ll,:)).^2).*weightvec'),2)/sum(weightvec);
            DF_denom_matrix{jj}(ll,:) = y;
            DF_matrix{jj}(ll,:) = DF_num_matrix{jj}(ll,:)./DF_denom_matrix{jj}(ll,:);
            DI_matrix{jj}(ll,:) = 10*log10(DF_matrix{jj}(ll,:));
        end
    end
%% When importing previously calculated data
elseif recalc == 2
    
    Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\results\ngaussvec.mat');
    ngaussvec = Dfile.ngaussvec;
    
    Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\results\DI_cell.mat');
    DI_matrix = Dfile.DI_matrix;
    
end
%% Figures

mm = length(ngaussvec)-1; %Place of DI_68 in DI_matrix

ref_DI = cell(length(frequency),1);
ref_DI_upper = cell(length(frequency),1);
ref_DI_lower = cell(length(frequency),1);
for hh = 1:length(frequency)
    ref_DI{hh} = (DI_matrix{hh}(mm)*ones(1,length(ngaussvec)));
    ref_DI_upper{hh} = (10*log10((10^((DI_matrix{hh}(mm))/10))*(10^(-3)+1)))*...
        ones(1,length(ngaussvec));
    ref_DI_lower{hh} = (10*log10((10^((DI_matrix{hh}(mm))/10))/(10^(-3)+1)))*...
        ones(1,length(ngaussvec));
end

figure(1)
fig1 = plot(ngaussvec,DI_matrix{1},'-o',ngaussvec,ref_DI{1},'-',...
    ngaussvec,ref_DI_upper{1},'--',ngaussvec,ref_DI_lower{1},'--');
set(fig1,'Linewidth',1)
fig1(1).Color = [0 0.4470 0.7410]; %blue
fig1(2).Color = [0.8500 0.3250 0.0980]; %red
fig1(3).Color = [0.8500 0.3250 0.0980]; %red
fig1(4).Color = [0.8500 0.3250 0.0980]; %red
g = legend('Simulated DI (10 kHz)','Reference DI_{68} (10 kHz)',...
    '0.0043 dB dev. from ref.');
set(g,'Location','best')
xlim([5 70])
% title({'Plot of simulated directivity index for different',...
%       'values of ngauss as input to EDtoolbox'})
xlabel({'Number of edge points, [1]'})
ylabel({'Directivity index, [dB]'})
set(gcf,'color','w')
grid on

figure(2)
fig2 = plot(ngaussvec,DI_matrix{2},'-o',ngaussvec,ref_DI{2},'-',...
    ngaussvec,ref_DI_upper{2},'--',ngaussvec,ref_DI_lower{2},'--');
set(fig2,'Linewidth',1)
fig2(1).Color = [0 0.4470 0.7410]; %blue
fig2(2).Color = [0.8500 0.3250 0.0980]; %red
fig2(3).Color = [0.8500 0.3250 0.0980]; %red
fig2(4).Color = [0.8500 0.3250 0.0980]; %red
g = legend('Simulated DI (20 kHz)','Reference DI_{68} (20 kHz)',...
    '0.0043 dB dev. from ref.');
set(g,'Location','best')
xlim([5 70])
% title({'Plot of simulated directivity index for different',...
%       'values of ngauss as input to EDtoolbox'})
xlabel({'Number of edge points, [1]'})
ylabel({'Directivity index, [dB]'})
set(gcf,'color','w')
grid on


