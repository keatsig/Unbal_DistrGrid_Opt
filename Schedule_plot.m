clear all; clc;

mins = repmat('00',12,1);
hr = num2str([0:2:22]');
x_time = strcat(hr,':',mins);

%% Data
case_name = 'case13_glm';
case_name = 'R1-12.47-3';
dat_load = xlsread(strcat('House_',case_name,'_profile.csv'))*1e3;
dat_PV = xlsread(strcat('PV_',case_name,'_profile.csv'))*1e3;
loc = strcat(case_name,'_DAY/min_PVUR/');
VU_LVUR = xlsread(strcat(loc,'LVUR_',case_name,'_soln.csv'));
VU_PVUR = xlsread(strcat(loc,'PVUR_',case_name,'_soln.csv'));
VU_VUF = xlsread(strcat(loc,'VUF_',case_name,'_soln.csv'));
Vm_soln = xlsread(strcat(loc,'Vm_',case_name,'_soln.csv'));

%% PV/Load/voltage profile plots
lw = 2; fs = 16;
figure('position', [10, 50, 1100, 950])
clf('reset')
subplot(3,1,1),
plot(dat_PV(:,1:96)','-o','linewidth',lw)
hold on; grid on;
plot(repmat([5],size(dat_PV,2),1),'r-.','linewidth',lw)
hold off;
xticks([1:8:96]); 
xticklabels(x_time);
xlabel('Time'); 
ylabel('PV Power [kW]');
xlim([1,96])
ylim([-0.2,5.5])
set(gca,'fontsize',fs)

subplot(3,1,2),
colors = distinguishable_colors(size(dat_load,1));
% colors = varycolor(size(dat_load,1));
grid on; hold on;
for i=1:size(dat_load,1)
    plot(dat_load(i,1:96),'-o','linewidth',lw,'color', colors(i,:))
end
hold off;
xticks([1:8:96]); 
xticklabels(x_time);
xlabel('Time'); 
ylabel('Load Power [kW]');
xlim([1,96])
ylim([-0.02,11])
set(gca,'fontsize',fs)

subplot(3,1,3),
colors = distinguishable_colors(size(Vm_soln,1));
hold on; grid on;
for i=1:size(Vm_soln,1)
    plot(Vm_soln(i,:),'-o','linewidth',lw,'color', colors(i,:))
end
plot(repmat([0.9],size(Vm_soln,2),1),'r-.','linewidth',lw)
plot(repmat([1.1],size(Vm_soln,2),1),'r-.','linewidth',lw)
hold off;
xticks([1:8:96]);  
xticklabels(x_time);
xlabel('Time'); 
ylabel('Voltage [pu]');
xlim([1,96])
ylim([min(0.95*min(min(Vm_soln)),0.85),max(1.05*max(max(Vm_soln)),1.15)])
set(gca,'fontsize',fs)

%% VU plots
[nb, nsch] = size(VU_LVUR);
for i=1:nsch
    VU_VUF(nb+1,i) = mean(VU_VUF(1:nb,i));
    VU_LVUR(nb+1,i) = mean(VU_LVUR(1:nb,i));
    VU_PVUR(nb+1,i) = mean(VU_PVUR(1:nb,i));
     VU_VUF(nb+2,i) = max(VU_VUF(1:nb,i));
    VU_LVUR(nb+2,i) = max(VU_LVUR(1:nb,i));
    VU_PVUR(nb+2,i) = max(VU_PVUR(1:nb,i));
end

lw = 2; lb = -0.2; ub = 2.5;
figure('position', [10, 50, 1100, 950])
subplot(3,1,1),
colors = distinguishable_colors(size(VU_VUF,1));
hold on; grid on;
for i=1:size(VU_VUF,1)
    plot(VU_VUF(i,:),'-o','linewidth',lw,'color', colors(i,:))
end 
hold off;
xticks([1:8:96]); 
xticklabels(x_time);
xlabel('Time'); 
ylabel('VUF [%]');
xlim([1,96])
ylim([lb,ub])
% title('VUF profile')
set(gca,'fontsize',fs)

subplot(3,1,2),
colors = distinguishable_colors(size(VU_LVUR,1));
hold on; grid on;
for i=1:size(VU_LVUR,1)
    plot(VU_LVUR(i,:),'-o','linewidth',lw,'color', colors(i,:))
end 
hold off; 
xticks([1:8:96]); 
xticklabels(x_time);
xlabel('Time'); 
ylabel('LVUR [%]');
xlim([1,96])
% ylim([lb,ub])
% title('LVUR profile')
set(gca,'fontsize',fs)

subplot(3,1,3),
colors = distinguishable_colors(size(VU_PVUR,1));
% colors = varycolor(size(VU_PVUR,1));
hold on; grid on;
for i=1:size(VU_PVUR,1)
    plot(VU_PVUR(i,:),'-o','linewidth',lw,'color', colors(i,:))
end 
hold off; 
xticks([1:8:96]); 
xticklabels(x_time);
xlabel('Time'); 
ylabel('PVUR [%]');
xlim([1,96])
ylim([lb,3])
% title('PVUR profile')
set(gca,'fontsize',fs)


