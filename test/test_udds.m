addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%  UDDS TEST  -  validate ecm_model.m on a driving-cycle dataset
%  1) Load UDDS.mat and CLEAN it (remove duplicate / non-increasing samples)
%  2) Build the SOC trajectory (coulomb counting, same eta/Q rule as the model)
%  3) Call ecm_model(t, I, SOC)  ->  model terminal voltage
%  4) Compare with the measured voltage: error and RMSE (in volts), plots

% UDDS is the Urban Dynamometer Driving Schedule, a standard city driving profile. 
% The cell is run with a current pattern that mimics driving: accelerations draw current (discharge), 
% braking or regeneration gives short charge bursts, and there are stops in between.

% The datset has a sampling of 0.1s and of 1s during the final rest
% The final rest contain multiple duplicate values
% The initial SOC during discharge is unknown 

clear ecm_model;            % reload the parameter tables (cached inside the model)

udds_file = 'UDDS.mat';
dt_min    = 0.01;        %samples closer than this to the previous one = duplicate
t_end_rest  = 60;        % last part of the file averaged as "rest voltage"
opts = struct('use_hyst', true, 'correct_H', true, 'rest_thr', 0.05); 

%  STEP 1: LOAD AND CLEAN THE DATA
%  It removes NaN files and the data where time stamps are closer than
%  0.01s(dt_min) are also removed so the time axis become strictly increasing.

u = load(udds_file);
t = u.Data.Times(:);   
V = u.Data.VoltageV(:);   
I = u.Data.CurrentA(:);

good = isfinite(t) & isfinite(V) & isfinite(I);          
t = t(good);  V = V(good);  I = I(good);

keep = [true; diff(t) > dt_min];                         % drop near-duplicate time stamps
t = t(keep);  V = V(keep);  I = I(keep);

%  SOC ESTIMATION
%  The UDDS file starts in the middle of a discharge (no rest at the start), so
%  the starting SOC cannot be read from the first voltage. 
%  There is a 90min rest at end of file so there is a proper rest for
%  voltage to drop closer to OCV (OCV - H).So for last 60s we take the average of
%  these volatges and find the SOC here via the look-up table and then
%  count back all the way to start to find the starting SOC..

o = load('OCV_SOC_table.mat', 'Q', 'eta');
Q = o.Q;   
eta = o.eta;

I_eff = I;
I_eff(I > 0) = eta * I(I > 0);
dSOC  = cumtrapz(t, I_eff) / 3600 / Q;                   % SOC change since the start

zg    = linspace(0, 1, 1001).';                  % SOC axis
Vrest = ecm_model((0:1000).', zeros(1001,1), zg, opts);   % rest V after discharge
        
[Vu, iu] = unique(Vrest);
Vend     = mean(V(t >= t(end) - t_end_rest));
SOC_end  = interp1(Vu, zg(iu), Vend, 'linear', 'extrap');
SOC0     = SOC_end - dSOC(end);
        
SOC = SOC0 + dSOC;
fprintf('SOC range during the test: %.1f %% ... %.1f %%\n', 100*min(SOC), 100*max(SOC));

% RUN THE MODEL
[Vm, parts] = ecm_model(t, I, SOC, opts);


% Here will validate the model prediction on the actual voltage data from
% the test file and will compute the RMSE for different intervals.
err   = Vm - V;
drive = abs(I) > opts.rest_thr;
rest  = ~drive & t > t(find(drive, 1, 'last'));          % final rest only

rmse = @(e) sqrt(mean(e.^2));
fprintf('\nRMSE (model - measured):\n');
fprintf('  Whole file : %.4f V\n', rmse(err));
fprintf('  Driving    : %.4f V\n', rmse(err(drive)));
if any(rest)
    fprintf('  Final rest : %.4f V\n', rmse(err(rest)));
end

% PLOTS
th = t / 3600;
figure('Name','UDDS: measured vs model','Color','w');
subplot(3,1,1);
plot(th, V, 'b', 'LineWidth', 1); hold on;
plot(th, Vm, 'r', 'LineWidth', 1); grid on;
ylabel('Voltage (V)'); legend('Measured','Model','Location','best');
title('UDDS test: terminal voltage');
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