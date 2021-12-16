% Script to calculate the directivity index for the simulations done on the
% box, including the original results, as well as the results for a maximum
% deviation in the sound velocity, the angle of transposition for the
% normal incidence, and the dimensions for the box.

% Written 14.11.21 by Isabel Berg
% Modified 30.11.21 to test for the impact of using appropriate kl axis for
% the different simulations
% Modified 01.12.21 to use "best index" method when calculating
% differences, and average kl value when plotting them.
% 
clear all
close all

frequency = (100:10:20000);
l = ((329.44 + 329.03 + 390.73)/3)/1000;
%l = 349.7341667/1000;
c = 344.7130;
c2 = 344.2873;
l4 = (((329.44+0.16) + (329.03+0.38) + (390.73+0.14))/3)/1000;
kl = (2*pi*frequency*l)/c;
kl2 = (2*pi*frequency*l)/c2;
kl4 = (2*pi*frequency*l4)/c;

%% Opening files

Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\Data\Box_simulation\box_sim_original.mat');
P_matrix_original = Dfile.tftot;

Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\Data\Box_simulation\box_sim_cdev.mat');
P_matrix_cdev = Dfile.tftot;

Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\Data\Box_simulation\box_sim_phivecdev.mat');
P_matrix_phivecdev = Dfile.tftot;

Dfile = load(...
    'C:\Users\isabe\Documents\MATLAB\Master\Data\Box_simulation\box_sim_adev.mat');
P_matrix_adev = Dfile.tftot;

P_matrix = cell(4,1);
P_matrix{1} = P_matrix_original;
P_matrix{2} = P_matrix_cdev;
P_matrix{3} = P_matrix_phivecdev;
P_matrix{4} = P_matrix_adev;

%% Calculations
weightvec = ones(size(P_matrix_original,2),1); %Weighting since only sim. 1:pi
weigthvec(1) = 0.5;
weightvec(end) = 0.5;
    
DF_num_matrix = cell(4,1);
DF_denom_matrix = cell(4,1);
DF_matrix = cell(4,1);
DI_matrix = cell(4,1);
for jj = 1:4
    for ll = 1:length(frequency)
        DF_num_matrix{jj}(ll) = abs((P_matrix{jj}(ll,1)).^2);
        y = sum((abs((P_matrix{jj}(ll,:)).^2).*weightvec'),2)/sum(weightvec);
        DF_denom_matrix{jj}(ll) = y;
        DF_matrix{jj}(ll) = DF_num_matrix{jj}(ll)./DF_denom_matrix{jj}(ll);
        DI_matrix{jj}(ll) = 10*log10(DF_matrix{jj}(ll));
    end
end


diff1 = DI_matrix{2}-DI_matrix{1}; % kl2, kl
diff2 = DI_matrix{3}-DI_matrix{1}; % kl, kl
diff3 = DI_matrix{4}-DI_matrix{1}; % kl4, kl

%% Best index:

indexnumber_in_klvec2 = findbestkavecindices(kl2,kl);
diff1_best_index = zeros(1,length(DI_matrix{1}));
for m = 1:length(DI_matrix{1})
    diff = DI_matrix{2}(indexnumber_in_klvec2(m)) - DI_matrix{1}(m);
    diff1_best_index(m) = diff;
end

indexnumber_in_klvec4 = findbestkavecindices(kl4,kl);
diff3_best_index = zeros(1,length(DI_matrix{1}));
for m = 1:length(DI_matrix{1})
    diff = DI_matrix{4}(indexnumber_in_klvec4(m)) - DI_matrix{1}(m);
    diff3_best_index(m) = diff;
end

%% Average kl for plotting difference: (code from Box_DI)
kl_avg = cell(1,3);
%Averaging over four series for each setup:
for l = 1:3
    if l == 1
        kl_avg{l} = (kl2 + kl)/2;
    elseif l == 3
        kl_avg{l} = (kl4 + kl)/2;
    else
        kl_avg{l} = kl;
    end
end

kl_diff_1 = kl2-kl;
kl_diff_3 = kl4-kl;

%% Figures

figure(1)
fig4 = semilogx(kl,DI_matrix{2},'-o',...
    kl,DI_matrix{3},'-o',kl,DI_matrix{4},'-o',kl,DI_matrix{1},'-o',...
    [kl(1466) kl(1466) kl(1466) kl(1466)],...
    [DI_matrix{2}(1466) DI_matrix{3}(1466) DI_matrix{4}(1466) DI_matrix{1}(1466)],...
    '*k');
set(fig4,'Linewidth',1)
fig4(1).Color = [0 0.4470 0.7410]; %blue
grid on
% title({'','Simulated DI using avg. temperature for',...
%     'four different measured series',''})
xlabel({'kl, [1]'})
xlim([11.96 12.01])
ylim([-0.07 0.05])
%xlim([94 94.16])
%ylim([4.495 4.52])
ylabel({'Directivity index, [dB]'})
labels = {'Max. deviation in c','0.2° off normal incidence',...
    'Max. dev. in box dimensions','Original simulation'};
neworder = [4 1 2 3]; % Re-order Legend
g = legend(fig4(neworder), labels(neworder));
set(g,'Location','best')
set(gcf,'color','w')

figure(2)
fig2 = semilogx(kl2,DI_matrix{2},'-o',...
    kl,DI_matrix{3},'-o',kl4,DI_matrix{4},'-o',kl,DI_matrix{1},'-o',...
    [kl2(1466) kl(1466) kl4(1466) kl(1466)],...
    [DI_matrix{2}(1466) DI_matrix{3}(1466) DI_matrix{4}(1466) DI_matrix{1}(1466)],...
    '*k');
set(fig2,'Linewidth',1)
fig2(1).Color = [0 0.4470 0.7410]; %blue
grid on
% title({'','Simulated DI using avg. temperature for',...
%     'four different measured series',''})
xlabel({'kl, [1]'})
xlim([11.96 12.01])
ylim([-0.07 0.05])
%xlim([94 94.16])
%ylim([4.495 4.52])
ylabel({'Directivity index, [dB]'})
labels = {'Max. deviation in c','0.2° off normal incidence',...
    'Max. dev. in box dimensions','Original simulation'};
neworder = [4 1 2 3]; % Re-order Legend
g = legend(fig2(neworder), labels(neworder));
set(g,'Location','best')
set(gcf,'color','w')

figure(3)
fig3 = semilogx(kl2,DI_matrix{2},'-o',...
    kl,DI_matrix{3},'-o',kl4,DI_matrix{4},'-o',kl,DI_matrix{1},'-o',...
    [kl2(1466) kl(1466) kl4(1466) kl(1466)],...
    [DI_matrix{2}(1466) DI_matrix{3}(1466) DI_matrix{4}(1466) DI_matrix{1}(1466)],...
    '*k');
set(fig3,'Linewidth',1)
fig3(1).Color = [0 0.4470 0.7410]; %blue
grid on
% title({'','Simulated DI using avg. temperature for',...
%     'four different measured series',''})
xlabel({'kl, [1]'})
% xlim([94 94.16])
% ylim([4.495 4.52])
ylabel({'Directivity index, [dB]'})
labels = {'Max. deviation in c','0.2° off normal incidence',...
    'Max. dev. in box dimensions','Original simulation'};
neworder = [4 1 2 3]; % Re-order Legend
g = legend(fig3(neworder), labels(neworder));
set(g,'Location','best')
set(gcf,'color','w')

figure(4)
fig6 = semilogx(kl_avg{1},diff1,'-',kl_avg{2},diff2,'-',kl_avg{3},diff3,'-');
set(fig6,'Linewidth',1)
%fig6(1).Color = [0 0.4470 0.7410]; %blue
grid on
xlim([0.6 135])
% title({'','Simulated DI using avg. temperature for',...
%     'four different measured series',''})
xlabel({'kl, [1]'})
ylabel({'Difference in directivity index, [dB]'})
g = legend('Max. deviation in c', '0.2° off normal incidence',...
    'Max. dev. in box dimensions');
set(g,'Location','best')
set(gcf,'color','w')


figure(5)
fig7 = semilogx(kl_avg{1},(diff1_best_index),'-',kl_avg{2},...
    (diff2),'-',kl_avg{3},(diff3_best_index),'-');
set(fig7,'Linewidth',1)
%fig6(1).Color = [0 0.4470 0.7410]; %blue
grid on
xlim([0.6 135])
% title({'','Simulated DI using avg. temperature for',...
%     'four different measured series',''})
xlabel({'kl, [1]'})
ylabel({'Difference in directivity index, [dB]'})
g = legend('Max. deviation in c', '0.2° off normal incidence',...
    'Max. dev. in box dimensions');
set(g,'Location','best')
set(gcf,'color','w')

% figure(6)
% fig8 = semilogx(kl,(diff1_best_index),'-',kl,...
%     (diff2),'-',kl,(diff3_best_index),'-');
% set(fig8,'Linewidth',1)
% %fig6(1).Color = [0 0.4470 0.7410]; %blue
% grid on
% xlim([0.6 135])
% % title({'','Simulated DI using avg. temperature for',...
% %     'four different measured series',''})
% xlabel({'kl, [1]'})
% ylabel({'Difference in directivity index, [dB]'})
% g = legend('Max. deviation in c', '0.2° off normal incidence',...
%     'Max. dev. in box dimensions');
% set(g,'Location','best')
% set(gcf,'color','w')

figure(7)
plot(kl_diff_1,'-o')
hold on
plot(kl_diff_3,'-o')
grid on
xlabel({'Sample no., [1]'})
ylabel({'Difference in kl number, [1]'})
g = legend('From dev. in c', 'From dev. in box dimensions');
set(g,'Location','best')
set(gcf,'color','w')
hold off
