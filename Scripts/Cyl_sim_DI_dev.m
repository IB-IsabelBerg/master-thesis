%Simulation of directivity index for four different sound speeds, using
%series 1-4, with the average of the temperature from each.

%Written 07.06.21 by Isabel Berg
%Edited 17.10.21 to add values for the radius as measured on the
%constructed cylinder object used in the anechoic chamber.
%Modified 01.12.21 to use average ka value when plotting differences.

clear all;
close all;


%%  Import Sound velocity file
%Opening data file for sound speed data
Dfile_1 = load('C:\Users\isabe\Documents\MATLAB\Master\Data\C_data_cyl.mat');
C_data = Dfile_1.C_data;

%Averaging for each measurement series
C_matrix_avg = zeros(8,1);
for k = 1:8
    C_matrix_avg(k,1) = mean(C_data(k,:),2);
end

%Average for all measurements on cylinder
C_avg = mean(C_matrix_avg);

%% DI by use of parameters from measurements on cylinder

%Frequency values, linearly distanced, arbitrarily chosen:
f =(50:10:20000);                     %Does this need a better resolution? Edit this after having looked at the SNR plots.
a = (0.20043)/2;%Cylinder radius (measured)   %pm 0.00012 for diameter 
M = 50; %Number of elements in the sum for the pressure equation
n_phi = 46; %Number of discrete angles for a semicircle around the cylinder
phi_vec = linspace(0,pi,n_phi).'; %Azimuth angle
%Weighting for only calculating half the values, but wanting to use a whole
%cycle around the cylinder:
weightvec = ones(length(phi_vec),1);
weigthvec(1) = 0.5;
weightvec(end) = 0.5;

lambda = C_avg./f; %Wave length of incoming plane wave
wavenum = (2*pi)./lambda; %Wave number
ka = wavenum.*a;

%DF_band_T1 = zeros (1,length(fcenter));
for i = 1:length(f)
    [P_a] = CylinderPressurev5(ka(i),M,phi_vec);
    if i == 1
        P_matrix = zeros(length(P_a),length(f));
    end
    P_matrix(:,i) = P_a;
end

DF_num = abs((P_matrix(n_phi,:)).^2);
DF_denom = zeros(1,length(DF_num));
for i = 1:length(f)
    y = sum((abs((P_matrix(:,i)).^2).*weightvec),1)/sum(weightvec);
    DF_denom(1,i) = y;
end
DF = DF_num./DF_denom;
DF_val = DF;

DI = 10*log10(DF);

% figure(1)
% fig1 = semilogx(ka,DI,'-');
% grid on
% %title({'Simulated Directivity Index for cylinder',' '})
% xlabel({'ka, [1]'})
% ylabel({'Directivity index, [dB]'})
% set(fig1,'LineWidth',1);
% g = legend('Simulated DI for cylinder');
% set(g,'Location','best');
% set(gcf,'color','w')

%% Sound velocity variation
%Comparing with the max value of C
C_max = max(max(C_data));
lambda_C2 = C_max./f; %Wave length of incoming plane wave
wavenum_C2 = (2*pi)./lambda_C2; %Wave number
ka_C2 = wavenum_C2.*a;

for i = 1:length(f)
    [P_a] = CylinderPressurev5(ka_C2(i),M,phi_vec);
    if i == 1
        P_matrix_C2 = zeros(length(P_a),length(f));
    end
    P_matrix_C2(:,i) = P_a;
end

DF_num = abs((P_matrix_C2(n_phi,:)).^2);
DF_denom = zeros(1,length(DF_num));
for i = 1:length(f)
    y = sum((abs((P_matrix_C2(:,i)).^2).*weightvec),1)/sum(weightvec);
    DF_denom(1,i) = y;
end
DF_C2 = DF_num./DF_denom;
DI_C2 = 10*log10(DF_C2);


%Finding difference between the two directivity indices:
deltadii_T = DI_C2 - DI;

%Testing with best index as well:
indexnumber_in_kavec1 = findbestkavecindices(ka_C2,ka);
deltadii_T_best_index = zeros(1,length(DI));
for m = 1:length(DI)
    diff = DI_C2(indexnumber_in_kavec1(m)) - DI(m);
    deltadii_T_best_index(m) = diff;
end

% figure(2)
% fig2 = semilogx(ka_C2,DI_C2,'-o',ka,DI,'-o');
% grid on
% title({'Simulated Directivity Index for cylinder with max C value',' '})
% xlabel({'ka, [1]'})
% ylabel({'Directivity index, [dB]'})
% set(fig2,'LineWidth',1);
% g = legend(['Max C = ', num2str(C_max),' m/s'],['Average C = ', ...
%     num2str(C_avg),' m/s']);
% set(g,'Location','best');
% set(gcf,'color','w')

% % % % figure(3)
% % % % %Directivity index difference of the two sound velocities:
% % % % fig3 = semilogx(ka,(deltadii_T),'-o');
% % % % %set(fig3,'LineWidth',1);
% % % % grid on
% % % % title({'Significance of sound velocity difference',''})
% % % % g = legend('Difference in DI for two sound velocities');
% % % % set(g,'Location','best');
% % % % xlabel({'ka, [1]'})
% % % % ylabel({'Directivity index, [dB]'})
% % % % set(gcf,'color','w')

% figure(4)
% %Directivity index difference of the two temperatures with best index:
% fig4 = semilogx(ka,(deltadii_T_best_index),'-o');
% %set(fig4,'LineWidth',1);
% grid on
% title({'Significance of sound velocity difference (best index)',''})
% g = legend('Difference in DI for two sound velocities (best index)');
% set(g,'Location','best');
% xlabel({'ka, [1]'})
% ylabel({'Directivity index, [dB]'})
% set(gcf,'color','w')

%% Slight variation in incident angle

% n_phi_NI = 46; %Number of discrete angles for a semicircle around the cylinder
phi_vec_NI = linspace(0.2,pi,n_phi).'; %Azimuth angle

for i = 1:length(f)
    [P_a] = CylinderPressurev5(ka(i),M,phi_vec_NI);
    if i == 1
        P_matrix_NI = zeros(length(P_a),length(f));
    end
    P_matrix_NI(:,i) = P_a;
end

%For 0.2 degrees towards the side: %Update angle value for cylinder        !!!
DF_num_NI = abs((P_matrix_NI((n_phi),:)).^2);
DF_NI = DF_num_NI./DF_denom;
DI_NI = 10*log10(DF_NI);

deltadii_NI = (DI_NI - DI);

% %Plot of 1 degree difference in incident angle
% figure(5)
% fig5 = semilogx(ka,(deltadii_NI),'-o');
% %set(fig5,'LineWidth',1);
% grid on
% title({'Significance of measured incident angle difference',''})
% xlabel({'ka, [1]'})
% ylabel({'Directivity index, [dB]'})
% g = legend('Difference between 0° incidence and 0.2° incidence');
% set(g,'Location','best');
% set(gcf,'color','w')

 %% Uncertainty in radius of cylinder 

a2 = (0.20043 + 0.00012)/2; %Simulating only radius larger than expected.
ka_a2 = wavenum.*a2;

for i = 1:length(f)
    [P_a] = CylinderPressurev5(ka_a2(i),M,phi_vec);
    if i == 1
        P_matrix_a2 = zeros(length(P_a),length(f));
    end
    P_matrix_a2(:,i) = P_a;
end

DF_nom_a2 = abs((P_matrix_a2(n_phi,:)).^2);
DF_denom_a2 = zeros(1,length(DF_nom_a2));
for i = 1:length(f)
    y = sum((abs((P_matrix_a2(:,i)).^2).*weightvec),1)/sum(weightvec);
    DF_denom_a2(1,i) = y;
end
DF_a2 = DF_nom_a2./DF_denom_a2;
DI_a2 = 10*log10(DF_a2);

%Finding difference between the two directivity indices:
deltadii_a = DI_a2 - DI;

%Using best index for comparison:
indexnumber_in_kavec1 = findbestkavecindices(ka_a2,ka);
deltadii_a_best_index = zeros(1,length(DI));
for m = 1:length(DI)
    diff = DI_a2(indexnumber_in_kavec1(m)) - DI(m);
    deltadii_a_best_index(m) = diff;
end


% % % % %Plot of measurement uncertainty difference in radius of infinite cylinder:
% % % % figure(6)
% % % % fig6 = semilogx(ka,(deltadii_a),'-o');
% % % % %set(fig6,'LineWidth',1);
% % % % grid on
% % % % title({'Significance of radius difference',''})
% % % % xlabel({'ka, [1]'})
% % % % ylabel({'Directivity index, [dB]'})
% % % % g = legend('Diff. measured radius and maximum meas. uncertainty');
% % % % set(g,'Location','best');
% % % % set(gcf,'color','w')

% figure(7)
% %Directivity index difference of the two temperatures with best index:
% fig7 = semilogx(ka,(deltadii_a_best_index),'-o');
% set(fig7,'LineWidth',1);
% grid on
% title({'Significance of radius difference (best index)',''})
% g = legend('Diff. measured radius and maximum meas. uncertainty (best index)');
% set(g,'Location','best');
% xlabel({'ka, [1]'})
% ylabel({'Directivity index, [dB]'})
% set(gcf,'color','w')

%% Final figures

ka_avg_T = (ka + ka_C2)/2;
ka_avg_a = (ka + ka_a2)/2;

figure(8)
fig8 = semilogx(ka_avg_T,(deltadii_T_best_index),'-',...
    ka,(deltadii_NI),'-',...
    ka_avg_a,(deltadii_a_best_index),'-');
set(fig8,'LineWidth',1)
grid on
g = legend('Max. deviation in c','0.2° off normal incidence',...
    'Max. dev. in radius');
set(g,'Location','best');
xlim([0.09 40])
ylim([-3.5 3]*(10^-3))
xlabel({'ka, [1]'})
ylabel({'Difference in directivity index, [dB]'})
set(gcf,'color','w')

figure(9)
fig9 = semilogx(ka_C2,DI_C2,'-',ka,DI_NI,'-',ka_a2,DI_a2,'-',ka,DI,'-');
set(fig9,'LineWidth',1)
grid on
%title({'Simulated Directivity Index for cylinder',' '})
xlim([0.09 40])
xlabel({'ka, [1]'})
ylabel({'Directivity index, [dB]'})
labels = {'Max. deviation in c',...
    '0.2° off normal incidence','Max. dev. in radius','Original simulation'};
neworder = [4 1 2 3]; % Re-order Legend
g = legend(fig9(neworder), labels(neworder));
set(g,'Location','best');
set(gcf,'color','w')



