% vG fits soil water retention model to measured data based on original van Genuchten
% IvG fits soil water retention model to measured data based on improved van Genuchten
% Suppose that alpha = vGx(1), theta_r = vGx(2), n = vGx(3), theta_s = vGx(4)
clear; clc; format longG;

th = readtable('C:\Users\ML\OneDrive - USU\Documents\MATLAB\Data_IvG.csv'); % Path of data
ts = 'Silt_Loam_UNSODA_3090';
rows_th = matches(th.Soil_sample,ts);
th1 = th(rows_th,:);
theta_v = th1.theta; 
h = th1.h; 

% van Genuchten model
vG = @(vGx,h) (1./(1+(vGx(1)*h).^(vGx(3)))).^(1-1/vGx(3))*(vGx(4)-vGx(2))+vGx(2); 

% Improved van Genuchten model
IvG = @(vGx,h) vGx(4)*(1./(1+(vGx(1)*h).^(vGx(3)))).^(1-1/vGx(3))...
    +vGx(2)*(log(10^6.8./h)).*(1-(1./(1+(vGx(1)*h).^(vGx(3)))).^(1-1/vGx(3)));

% Define initial parameter values and upper and lower bounds
vGx0 = [0.01; 0.001; 1.5; 0.01];
lb_vG = [0; 0; 1; 0];
% ub_vG = [1; 0.45/20; 10; 1]; % For Clay soil
ub_vG = [1; 1; 10; 1];

% Fit model to measured data using lsqcurvefit function
[vGx] = lsqcurvefit(vG,vGx0,h,theta_v,lb_vG,ub_vG)
[IvGx] = lsqcurvefit(IvG,vGx0,h,theta_v,lb_vG,ub_vG)

% Calculate root mean squared error (RMSE)
RMSE = {'SWRC Model', 'RMSE', 'Unit';...    % Build cell array
        'van Genuchten', sqrt(mean((vG(vGx,h)-theta_v).^2)),'cm3/cm3';...
        'Improved van Genuchten', sqrt(mean((IvG(IvGx,h)-theta_v).^2)),'cm3/cm3'}

% R2_vG = rsquare(vG(vGx,h),theta_v)
% R2_IvG = rsquare(IvG(IvGx,h),theta_v)

% %%% ----------    Create Chart with Two y-Axes  ---------------------------
h_min = min(h)/10; h_max = 10^6.8; % upper and lower limits for plotting
h_plot = h_min:1:h_max; 

figure(1);
set(figure(1),'defaulttextinterpreter','latex');
plot(log10(h_plot),vG(vGx,h_plot),'r-.','LineWidth',2.5); hold on % VG
plot(log10(h_plot),IvG(IvGx,h_plot),'b-','LineWidth',2); hold on % IVG
scatter(log10(h),theta_v,'ko','SizeData',80);

xlabel('- Pressure head, \bf{$h~(\tt{cm})$}','FontSize',14);
ylabel('Volumetric water content, \bf{$\theta~(\tt{cm^{3}}~\tt{cm^{-3}}) $}','FontSize',14);
legend('van Genuchten','Improved van Genuchten','Measured values')
text(0.2,0.04,'Silt Loam (#3090)','color','k','FontSize',14)
% xticklabels({'10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}'})
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}'})
% xticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}','10^{6}','10^{8}'})
