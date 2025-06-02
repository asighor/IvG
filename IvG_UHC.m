% Prediction of the unsaturated hydraulic conductivity based on vG and IvG

clear; clc; format long;

% Load data from Excel
filename = 'Your data path\adelanto_loam_IvG_UHC.xls';

Data = readmatrix(filename);

ks = 5.0; % Adelanto Loam
% ks = 12.0; % Pachappa Loam
% ks = 673.0; % Shonai Sand
% ks = 3.5; % silty clay Canning

% Extract relevant columns
h = Data(1:20,[1]); theta = Data(1:20,[2]); th = Data(1:6,[3]); Km = Data(1:6,[4]); % Adelanto Loam
% h = Data(1:23,[1]); theta = Data(1:23,[2]); th = Data(1:10,[3]); Km = Data(1:10,[4]); % Pachappa Loam
% h = Data(1:31,[1]); theta = Data(1:31,[2]); th = Data(1:67,[3]); Km = Data(1:67,[4]); % Shonai sand
% h = Data(1:10,[1]); theta = Data(1:10,[2]); th = Data(1:11,[3]); Km = Data(1:11,[4]); % silty clay Canning

% Define van Genuchten and improved van Genuchten functions
fun_VG = @(params,h) params(2) + (params(4)-params(2)) .* (1 ./ (1 + (params(1)*h).^params(3))).^(1 - 1./params(3));

fun_IVG = @(params,h) params(4)*(1./(1+(params(1)*h).^params(3))).^(1 - 1/params(3)) + ...
    params(2)*log(10^6.8 ./ h).*(1 - (1./(1+(params(1)*h).^params(3))).^(1 - 1/params(3)));

% Initial guesses and bounds for vG
x0_vg = [0.01; 0.001; 1.5; 0.01];
lb_vg = [0; 0; 1; 0];
ub_vg = [1; 1; 10; 1];

% Fit vG model
vg_model = @(p,h) fun_VG(p,h);
vg_params = lsqcurvefit(vg_model, x0_vg, h, theta, lb_vg, ub_vg)

% Initial guesses and bounds for IvG
x0_ivg = [0.01; 0.001; 1.5; 0.01];
lb_ivg = [0; 0; 1; 0];
ub_ivg = [1, 1, 10, 1];

% Fit IvG model
ivg_model = @(p,h) fun_IVG(p,h);
ivg_params = lsqcurvefit(ivg_model, x0_ivg, h, theta, lb_ivg, ub_ivg);

% Assign fitted parameters for IvG
alpha = ivg_params(1); 
n = ivg_params(3); 
m = 1 - 1/n;
qr = ivg_params(2); 
qs = ivg_params(4); 

% Assign fitted parameters for vG
alpha1 = vg_params(1); 
n1 = vg_params(3); 
m1 = 1 - 1/n1;
qr1 = vg_params(2); 
qs1 = vg_params(4);

% Solve IvG numerically using Gaussian quadrature
N = 1;
[X,W] = lgwt(N, -1, 1);

th_IvG = linspace(qr,qs+0.01,100);

for i = 1:numel(th_IvG)
    S = abs((th_IvG(i) - qr)/(qs - qr));

    for j = 1:N
        S1 = (S/2)*X(j) + (S/2);
        ht_func = @(h) qs/(qs-qr)*(1/(1+(alpha*h)^n)^m) + qr/(qs-qr)*(log(10^6.8/h)*(1-(1/(1+(alpha*h)^n)^m))-1) - S1;
        htt_func = @(h) qs/(qs-qr)*(1/(1+(alpha*h)^n)^m) + qr/(qs-qr)*(log(10^6.8/h)*(1-(1/(1+(alpha*h)^n)^m))-1) - (1/2*X(j) + 1/2);

        % Ensure endpoints of interval bracket a root
        ht(j) = fzero(ht_func, fzero_interval(ht_func, 1e-6, 1e+99));
        htt(j) = fzero(htt_func, fzero_interval(htt_func, 1e-6, 1e+99));
    end

    ft = (S/2)*sum(W./ht);
    f1 = 1/2*sum(W./htt);
    K_IvG(i) = ks*S^0.5*(ft/f1)^2;
end

th_vG = linspace(qr1,qs1+0.01,100);

for ii = 1:numel(th_vG)
    S = abs((th_vG(ii) - qr1)/(qs1 - qr1));
    K_vG(ii) = real(ks*S^0.5*(1 - (1 - S^(1/m1))^m1)^2);
end

% Plotting
figure;
semilogy(th, Km, 'ko', 'MarkerFaceColor','k'); hold on;
semilogy(th_IvG, K_IvG, 'b-','LineWidth',1);
semilogy(th_vG, K_vG, 'r-','LineWidth',1);
legend('Data', 'IVG Fit', 'VG Fit','Location','northoutside','Orientation','horizontal');
xlabel('\bf{Volumetric water content, \theta (cm^{3} cm^{-3})}'); 
ylabel('\bf{Hydraulic conductivity, K (cm d^{-1})}'); 
grid on; box on;

function interval = fzero_interval(f, a, b)
    % Expand interval until it brackets a root
    max_iter = 1000;
    factor = 1.5;
    for k = 1:max_iter
        if sign(f(a)) ~= sign(f(b))
            interval = [a, b];
            return;
        end
        a = max(a / factor, 1e-12);
        b = min(b * factor, 1e+12);
    end
    error('No sign change found in interval. Consider adjusting bounds or checking function.');
end

function [x, w] = lgwt(N, a, b)
    beta = 0.5 ./ sqrt(1 - (2*(1:N-1)).^(-2));
    T = diag(beta,1) + diag(beta,-1);
    [V, D] = eig(T);
    x = diag(D);
    [x, i] = sort(x);
    w = 2 * (V(1,i)').^2;
    x = 0.5 * ((b - a) * x + (b + a));
    w = 0.5 * (b - a) * w;
end
