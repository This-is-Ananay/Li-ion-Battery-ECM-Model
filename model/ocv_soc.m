%  PSEUDO-OCV ANALYSIS
%  Dataset: % https://data.mendeley.com/datasets/ycx459r5c3/2
%  charge current >0 , discharge <0

chg_file = load('C20_CHARGE.mat');   
dis_file = load('C20_DISCHARGE.mat');

% Coloumn vector of Time, Voltage, Current
t_chg = chg_file.Data.Times(:);      
V_chg_unique = chg_file.Data.VoltageV(:);   
I_chg = chg_file.Data.CurrentA(:);   

t_dis = dis_file.Data.Times(:);
V_dis_unique = dis_file.Data.VoltageV(:);
I_dis = dis_file.Data.CurrentA(:);   

% Capacity
Q_charge = cumtrapz(t_chg, I_chg) / 3600;     
Q_discharge = cumtrapz(t_dis, I_dis) / 3600;     

Q_chg_total = Q_charge(end);               % total Ah during charge
Q_dis_total = abs(Q_discharge(end));       % total Ah during discharge


%  Capacity
%  The charge starts near 2.994 V but the discharge ends at 2.5 V,
%  so the discharge covers a slightly LARGER voltage window than the charge.
%  Dividing the two raw totals would give eta > 1 (not possible)
%  So we compare in same range between the top voltage (~4.19 V) and the
%  voltage where the charge start.

[V_dis_sorted, ia] = unique(V_dis_unique);              
Qdis_abs_sorted = abs(Q_discharge(ia));            
Q_dis_window = interp1(V_dis_sorted, Qdis_abs_sorted, V_chg_unique(1));  

eta = Q_dis_window / Q_chg_total;          % effective coulombic efficiency
Q = Q_dis_total;                           % cell capacity [Ah] (full discharge)

fprintf('Charge capacity    : %.4f Ah\n', Q_chg_total);
fprintf('Discharge (window) : %.4f Ah\n', Q_dis_window);
fprintf('Discharge capacity : %.4f Ah\n', Q_dis_total);
fprintf('Coulombic eff. eta : %.5f (%.3f %%)\n', eta, 100*eta);

% State Of Charge Estimation
% Starting SOC is assumed as 0 during Charging and 1 during Discharging

SOC_dis = 1 + (1 * Q_discharge)/ Q;           
SOC_chg = 0 + (eta * Q_charge)/ Q;         

%  We interpolate VOLTAGE as a function of SOC (not SOC as a function of
%  voltage), because SOC is always monotonic while voltage can be flat or
%  noisy (especially in the CV tail of the charge).

% We need to interpolate both the Charge and Discharge SOC relationship
% before averaging them out. interp1 requires the SOC to be strictly
% increasing and unique so in both of the sets will use unique() to remove
% the duplicate values and sort them in ascending order.

[SOC_chg_unique, idx_c] = unique(SOC_chg);     % idx_c is storing actual index of the sorted value from previous data.
V_chg_unique = V_chg_unique(idx_c);            % using this we get that actual voltage corresponding to the SOC

[SOC_dis_unique, idx_d] = unique(SOC_dis);
V_dis_unique = V_dis_unique(idx_d);

% INTERPOLATION + EXTRAPOLATION
% Will also extrapolate the SOC to 100% if they are ending before it but
% these values are less reliable as voltage is steep near the end.

SOC_grid = linspace(0, 1, 1001).';           % linspace gives row vector so .' transpose it to coloumn.
V_chg_i = interp1(SOC_chg_unique, V_chg_unique, SOC_grid, 'linear', 'extrap');
V_dis_i = interp1(SOC_dis_unique, V_dis_unique, SOC_grid, 'linear', 'extrap');

%  These relationships are with the terminal voltage as the datset is
%  pseudo and don't provide any rest case so will compute the OCV by
%  averaging  these.
%       OCV(SOC) = ( V_chg(SOC) + V_dis(SOC) ) / 2

OCV = (V_chg_i + V_dis_i) / 2;
half_drop_est = (V_chg_i - V_dis_i) / 2;       % half-gap = IR/RC drop + hysteresis


figure('Name','OCV-SOC','Color','w');
plot(SOC_grid*100, V_chg_i, 'r--', 'LineWidth', 1.0); hold on;
plot(SOC_grid*100, V_dis_i, 'b--', 'LineWidth', 1.0);
plot(SOC_grid*100, OCV,     'k',   'LineWidth', 2.0);
grid on;
xlabel('SOC (%)');
ylabel('Voltage (V)');
title('Pseudo-OCV: average of charge and discharge curves');
legend('Charge (interp.)','Discharge (interp.)','OCV (average)','Location','southeast');


% STATIC HYSTERESIS VOLTAGE APPROXIMATION
% Here we are computing only Static Hysteresis
% See half gap is not exactly hysteresis it is mixture of hysteresis and
% ohmic + polarisation drop.

%  H(SOC) = ( V_chg(SOC) - V_dis(SOC) ) / 2  -  I_rate * R_dc

%  How H is cleaned:
%    only the SOC range where BOTH curves really have data is used (no
%    extrapolated voltages, they are steep and unreliable near the ends);
%    outside that range H is held at its nearest valid value
%    light moving-average smoothing (the difference of two curves is noisy)
%    H is never allowed to be negative

% Currently we use R_dc=0 and so will assume whole half gap is the
% hysteresis drop and obtain the relationship later on in ecm_model.m when
% we have the Other parameters will subtaract the drop accoridngly and use
% the correct value for further prediction...
R_dc      = 0;        
I_rate    = Q/20;     % (not useful)
smooth_n  = 21;       

soc_lo = max(min(SOC_chg_unique), min(SOC_dis_unique));      % overlap of the two curves
soc_hi = min(max(SOC_chg_unique), max(SOC_dis_unique));
in_ov  = SOC_grid >= soc_lo & SOC_grid <= soc_hi;

Vc_overlap  = interp1(SOC_chg_unique, V_chg_unique, SOC_grid(in_ov), 'linear');   % no extrap
Vd_overlap  = interp1(SOC_dis_unique, V_dis_unique, SOC_grid(in_ov), 'linear');

H_overlap = (Vc_overlap - Vd_overlap) / 2 - I_rate * R_dc;        % half-gap minus IR/RC drop
H_overlap = movmean(H_overlap, smooth_n, 'Endpoints', 'shrink');
H_overlap = max(H_overlap, 0);

% Interpolate and extrapolate the Hysteresis voltage
H = interp1(SOC_grid(in_ov), H_overlap, SOC_grid, 'nearest', 'extrap');

V_chg_h = OCV + H;
V_dis_h = OCV - H;

fprintf('Hysteresis H: mean %.1f mV, max %.1f mV (at SOC = %.1f %%)\n', 1e3*mean(H), 1e3*max(H), 100*SOC_grid(find(H==max(H),1)));

%  HYSTERESIS PLOT
figure('Name','Hysteresis','Color','w');
x = SOC_grid*100;
fill([x; flipud(x)], [V_chg_h; flipud(V_dis_h)], [0.2 0.45 0.9],'FaceAlpha', 0.30, 'EdgeColor', 'none');
hold on;
plot(x, V_chg_h, 'r', 'LineWidth', 1.5);
plot(x, V_dis_h, 'b', 'LineWidth', 1.5);
plot(x, OCV, 'k--','LineWidth', 1.5);
grid on;
xlabel('SOC (%)'); ylabel('Voltage (V)');
title('Voltage hysteresis between charge and discharge');
legend('Hysteresis band (2H)','Charge branch (OCV + H)','Discharge branch (OCV - H)','OCV (average)','Location','southeast');

% SOC-HYSTERESIS LOOKUP TABLE
Hyst_table = table(SOC_grid, H, 'VariableNames', {'SOC','H_Voltage'});   

% SOC-OCV LOOKUP TABLE
OCV_table = table(SOC_grid, OCV, 'VariableNames', {'SOC','OCV_Voltage'});

% SAVE THEM AS OCV_SOC_TABLE FILE
save(fullfile(fileparts(which('setup_paths')),'lookup','OCV_SOC_table.mat'),'OCV_table', 'Hyst_table', 'eta', 'Q');



%fprintf('\n       OCV-SOC LOOKUP TABLE\n');
%disp(OCV_table);
%fprintf('\n       HYSTERESIS LOOKUP TABLE\n');
%disp(Hyst_table);

