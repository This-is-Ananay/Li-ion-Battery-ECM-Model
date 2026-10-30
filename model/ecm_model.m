% MODEL FUNCTION 
% This contains a function that needs Time , Current , SOC and opts
% (optional struct):
%           use_hyst  true/false to use hysteresis (default = True)
%           correct_H true/false to use the corrected value of hysteresis
%           rest_thr   below this keep last direction s  (s>1 for charge)

% It will give the terminal voltage V as the output

%  Model:  V = OCV(SOC) + s*H(SOC) + R0*I + R1*I_R1
%           I_R1[k+1] = a*I_R1[k] + (1-a)*I[k],   a = exp(-dt/(R1*C1)) ( in
%           discrete form)
%           s = sign(I), last non-zero value is held at rest
%           R0, R1, C1: charge table if s > 0, discharge table if s < 0

% DOUBLE COUNTING is fixed here and will remove the polarisation part form
% the hysteresis we computed in ocv_soc.m to get the corrected hystersis
% voltage and will use that for prediction. With correct_H = true the model
% removes it H_used = max( H_saved - (I_c20)*(R0+R1), 0 ),  with R0+R1 
% taken as the average of the charge and discharge tables at each SOC.

% Here the charge - dischareg current in the datset is 2.9001 so will use
% that value only, with a limitation of CV charging short period.

function [V, parts] = ecm_model(t, I, SOC, opts)

if nargin < 3
    help ecm_model
    V = [];  parts = [];
    return
end

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'use_hyst'),  opts.use_hyst  = true; end
if ~isfield(opts, 'correct_H'), opts.correct_H = true; end
if ~isfield(opts, 'rest_thr'),  opts.rest_thr  = 0.05; end

% ---- load parameters (cached between calls) ----
persistent M
if isempty(M)
    o = load('OCV_SOC_table.mat');           % OCV_table, Hyst_table, eta, Q
    h = load('HPPC_params.mat');             % LUT
    L = h.LUT;
    M.Q    = o.Q;          M.eta = o.eta;
    M.OCV  = o.OCV_table.OCV_Voltage;              % uniform 0..1 SOC grid, 1001 pts
    M.Hraw = o.Hyst_table.H_Voltage;
    M.R0c  = L.R0_charge_Ohm;     M.R0d = L.R0_discharge_Ohm;
    M.R1c  = L.R1_charge_Ohm;     M.R1d = L.R1_discharge_Ohm;
    M.C1c  = L.C1_charge_F;       M.C1d = L.C1_discharge_F;
    Rdc    = 0.5*((M.R0c + M.R1c) + (M.R0d + M.R1d));      % [Ohm] at each SOC
    I_c20 = 2.9001;
    M.Hcor = max(M.Hraw - (I_c20) * Rdc, 0);              % IR/RC drop removed
end
if opts.correct_H, H = M.Hcor; else, H = M.Hraw; end

t = t(:);  I = I(:);  n = numel(t);

soc_given = numel(SOC) == n && n > 1;
if soc_given
    SOC = SOC(:);
else
    s0 = SOC(1);  SOC = zeros(n, 1);  SOC(1) = s0;     % tracked below
end

V    = zeros(n, 1);
I_R1 = zeros(n, 1);
pOCV = zeros(n, 1);  pH = zeros(n, 1);  pR0 = zeros(n, 1);  pRC = zeros(n, 1);
sv   = zeros(n, 1);
s = -1;                                      % last non-zero current direction

for k = 1:n
    z = min(max(SOC(k), 0), 1);
    if abs(I(k)) > opts.rest_thr, s = sign(I(k)); end

    if s > 0
        R0 = lut(M.R0c, z);  R1 = lut(M.R1c, z);  C1 = lut(M.C1c, z);
    else
        R0 = lut(M.R0d, z);  R1 = lut(M.R1d, z);  C1 = lut(M.C1d, z);
    end

    pOCV(k) = lut(M.OCV, z);
    if opts.use_hyst, pH(k) = s * lut(H, z); end
    pR0(k)  = R0 * I(k);
    pRC(k)  = R1 * I_R1(k);
    sv(k)   = s;
    V(k)    = pOCV(k) + pH(k) + pR0(k) + pRC(k);

    if k < n
        dt = t(k+1) - t(k);
        a  = exp(-dt / (R1*C1));
        I_R1(k+1) = a*I_R1(k) + (1-a)*I(k);
        if ~soc_given
            if I(k) > 0, e = M.eta; else, e = 1; end
            SOC(k+1) = SOC(k) + e*I(k)*dt / (3600*M.Q);
        end
    end
end

parts = struct('OCV', pOCV, 'hyst', pH, 'R0', pR0, 'RC', pRC,'I_R1', I_R1, 'SOC', SOC, 's', sv);
end

% FUNCTION 
% If the given SOC to model is not exactly in the LOOK-UP table it
% INTERPOLATE the values of parameters for that SOC.
function y = lut(tab, z)
    p = z * (numel(tab) - 1) + 1;            % fractional index
    i = min(floor(p), numel(tab) - 1);
    y = tab(i) + (p - i) * (tab(i+1) - tab(i));
end