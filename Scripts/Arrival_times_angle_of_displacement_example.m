%Script to illustrate arrival times and the process used to find the angle
%of displacement for achieving normal incidence.

%Written 17.11.21 by Isabel Berg




Dfile = load('C:\Users\isabe\Documents\MATLAB\Master\Data\Angle of incidence\Arrival_times_18');
Arrivaltimes = Dfile.table;

Dfile = load('C:\Users\isabe\Documents\MATLAB\Master\Data\Angle of incidence\Arrival_times_13');
Arrivaltimes2 = Dfile.table;

ampl = Arrivaltimes2(:,7:12);

maxampl = max(ampl);
minampl = min(ampl);
maxmaxampl = max([ maxampl abs(minampl)]);

normArrivaltimes2 = Arrivaltimes2/maxmaxampl;


% figure(1)
% subplot(2,1,1)
% plot(Arrivaltimes(:,1),'-')
% grid on
% xlim([285 300])
% ylim([-245 210])
% hold on
% plot(Arrivaltimes(:,2),'-')
% g = legend('-90°','+90°');
% set(g,'Location','best')
% set(gcf,'color','w')
% hold off
% 
% subplot(2,1,2)
% plot(Arrivaltimes(:,3),'-')
% grid on
% xlim([285 300])
% ylim([-245 210])
% hold on
% plot(Arrivaltimes(:,4),'-')
% g = legend('-90°','+90°');
% set(g,'Location','best')
% set(gcf,'color','w')
% hold off

figure(1)
% subplot(4,1,1)
% fig1 = plot(Arrivaltimes2(:,1),'-');
% grid on
% xlim([285 300])
% ylim([-245 235])
% hold on
% fig2 = plot(Arrivaltimes2(:,2),'-'); 
% fig3 = plot(Arrivaltimes2(:,3),'-');
% fig4 = plot(Arrivaltimes2(:,4),'-');
% fig6 = plot(Arrivaltimes2(:,5),'-');
% fig5 = plot(Arrivaltimes2(:,6),'-');
% fig1.Color = [0 0 1]; %standard blue
% fig2.Color = [1 0 0]; %standard red
% fig3.Color = [0.3010 0.7450 0.9330]; %light blue
% fig4.Color = [0.6350 0.0780 0.1840]; %dark red
% fig6.Color = [0 0.4470 0.7410]; %blue
% fig5.Color = [0.8500 0.3250 0.0980]; %orange
% g = legend('-90° (-7.4°)','+90° (-7.4°)','-90° (-8.5°)','+90° (-8.5°)',...
%     '-90° (-10.0°)','+90° (-10.0°)');
% set(g,'Location','best')
% set(gcf,'color','w')
% hold off

subplot(3,1,1);
fig1 = plot(normArrivaltimes2(:,7),'-');
set(fig1,'Linewidth',1)
grid on
xlim([285 300])
%ylim([-245 235])
hold on
fig2 = plot(normArrivaltimes2(:,8),'-');
set(fig2,'Linewidth',1)
g = legend('-90° (-12.0°)','+90° (-12.0°)');
%set(g,'Location','best')
set(gcf,'color','w')
set(gca,'fontsize',10)
%ylabel({'Impulse response'})
hold off

subplot(3,1,2)
fig3 = plot(normArrivaltimes2(:,11),'-');
set(fig3,'Linewidth',1)
grid on
xlim([285 300])
%ylim([-245 235])
hold on
fig4 = plot(normArrivaltimes2(:,12),'-');
set(fig4,'Linewidth',1)
g = legend('-90° (-12.1°)','+90° (-12.1°)');
%set(g,'Location','best')
set(gcf,'color','w')
set(gca,'fontsize',10)
ylabel({'Impulse response (normalized), [1]'})
hold off

subplot(3,1,3)
fig5 = plot(normArrivaltimes2(:,9),'-');
set(fig5,'Linewidth',1)
grid on
xlim([285 300])
%ylim([-245 235])
hold on
fig6 = plot(normArrivaltimes2(:,10),'-');
set(fig6,'Linewidth',1)
g = legend('-90° (-12.2°)','+90° (-12.2°)');
%set(g,'Location','best')
set(gcf,'color','w')
set(gca,'fontsize',10)
xlabel({'Sample no., [1]'})
%ylabel({'Impulse response'})
hold off



% figure(2)
% plot((1:20000),Arrivaltimes(:,1),(1:20000),Arrivaltimes(:,2))
% grid on
% xlim([285 300])
% ylim([-245 210])
% g = legend('-90°','+90°');
% set(g,'Location','best')
% set(gcf,'color','w')


% filestart = 'M08';
% fileend = '_S01_R01.etx';
% nfiles = 12;
% nuseful = 20000;
% readorigfiles = 1;
% 
% if readorigfiles == 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % Read IR files for set 1
% 
% for ii = 1:nfiles
%        filename = [filestart,int2str(58+ii),fileend];
%        FID=fopen(filename);
%        cdata1=textscan(FID,'%f%f','delimiter',',', 'HeaderLines', 22 );
%        fclose(FID);    
%        ir = cdata1{2};
%        tvecstart = cdata1{1}(1:2);
% 	fs = 1/(tvecstart(2)-tvecstart(1));
% 
% 	ir = ir(1:nuseful);
% 
% 	if ii == 1
% 		allirs_set1 = zeros(length(ir),nfiles);
% 	end	
% 	  allirs_set1(:,ii) = ir;
% end
% end
% 
% %Used in process of saving data in .mat-files:
% table = allirs_set1(1:20000,:);
