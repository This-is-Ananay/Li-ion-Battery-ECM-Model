addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% US06 is a more aggresive driving test particulary done on Highways.
% It has discharge upto -60A , regenerative charge pulses upto +27A
% The ending is a rest period of 90min 
% The driving ends at 2.5V cutoff voltage and has a very low SOC at end
% So will compute RMSE for whole driving and for driving greter than 0.1
% soc as the low soc at end is major culprit for higher error.

clear ecm_model;            % reload the parameter tables (cached inside the model)

us06_file = 'US06.mat';
dt_min    = 0.01;        % samples closer than this to the previous one = duplicate
t_end_rest  = 60;        % last part of the file averaged as "rest voltage"
soc_low     = 0.10;      % SOC below this = "low-SOC end" of the drive (extra RMSE line)
opts = struct('use_hyst', true, 'correct_H', true, 'rest_thr', 0.05);  

% LOAD AND CLEAN
u = load(us06_file);
t = u.Data.Times(:);   V = u.Data.VoltageV(:);   I = u.Data.CurrentA(:);

good = isfinite(t) & isfinite(V) & isfinite(I);         
t = t(good);  V = V(good);  I = I(good);

keep = [true; diff(t) > dt_min];                     % drop near-duplicate time stamps
t = t(keep);  V = V(keep);  I = I(keep);


%  SOC ESTIMATION
o = load('OCV_SOC_table.mat', 'Q', 'eta');
Q = o.Q;   eta = o.eta; 

I_eff = I;
I_eff(I > 0) = eta * I(I > 0);
dSOC  = cumtrapz(t, I_eff) / 3600 / Q;                   % SOC change since the start

zg = linspace(0, 1, 1001).';                              % SOC axis
Vrest = ecm_model((0:1000).', zeros(1001,1), zg, opts);   % rest V after discharge
[Vu, iu] = unique(Vrest);
Vend    = mean(V(t >= t(end) - t_end_rest));
SOC_end  = interp1(Vu, zg(iu), Vend, 'linear', 'extrap');
SOC0    = SOC_end - dSOC(end);
       
SOC = SOC0 + dSOC;
fprintf('SOC range during the test: %.1f %% ... %.1f %%\n', 100*min(SOC), 100*max(SOC));

% MODEL
[Vm, parts] = ecm_model(t, I, SOC, opts);

%ERROR
err   = Vm - V;
drive = abs(I) > opts.rest_thr;
rest  = ~drive & t > t(find(drive, 1, 'last'));          % final rest only
drive_hi = drive & SOC > soc_low;                        % driving, away from the empty end
drive_lo = drive & SOC <= soc_low;

rmse = @(e) sqrt(mean(e.^2));
fprintf('\nRMSE (model - measured):\n');
fprintf('  Whole file          : %.4f V\n', rmse(err));
fprintf('  Driving             : %.4f V\n', rmse(err(drive)));
if any(drive_hi)
    fprintf('  Driving, SOC > %2.0f %% : %.4f V\n', 100*soc_low, rmse(err(drive_hi)));
end
if any(drive_lo)
    fprintf('  Driving, SOC <= %2.0f %%: %.4f V\n', 100*soc_low, rmse(err(drive_lo)));
end
if any(rest)
    fprintf('  Final rest          : %.4f V\n', rmse(err(rest)));
end

% PLOTS
th = t / 3600;
figure('Name','US06: measured vs model','Color','w');
subplot(3,1,1);
plot(th, V, 'b', 'LineWidth', 1); hold on;
plot(th, Vm, 'r', 'LineWidth', 1); grid on;
ylabel('Voltage (V)'); legend('Measured','Model','Location','best');
title('US06 test: terminal voltage');
subplot(3,1,2);
plot(th, err, 'k', 'LineWidth', 0.8); grid on;
ylabel('Error (V)'); title('Model - measured');
subplot(3,1,3);
plot(th, I, 'Color', [0.1 0.5 0.1]); grid on;
xlabel('Time (h)'); ylabel('Current (A)');

% Error vs SOC
figure('Name','Error vs SOC','Color','w');
plot(100*SOC, err, '.', 'MarkerSize', 3); grid on;
xlabel('SOC (%)'); ylabel('Error (V)');
title('Voltage error vs SOC (drive + rest)');