% HPPC TEST TO COMPUTE 1-RC MODEL PARAMETERS

% R0 : instant voltage jump when the current starts/stops
% RC : exponential fit of the voltage relaxation after the pulse (over 5RC)
% R1 : relaxation voltage dV1, corrected for the short pulse
% C1 : RC / R1

% See HPPC is a Pulse test with several short pulses and rest period
% occuring, so main task is to identify the pulse.
k_rc          = 5;     % taking a settled relaxation of 5RC
rest_thr      = 0.5;   
pulse_min_dur = 5;     
pulse_max_dur = 30;    

% Load Data
hppc = load('HPPC.mat');
t = hppc.Data.Times(:);        
V = hppc.Data.VoltageV(:);     
I = hppc.Data.CurrentA(:);     

ocv = load('OCV_SOC_table.mat', 'Q');
Q   = ocv.Q;                   % Capacity (full discharge)

% Computing SOC as dataset don't have any related coloumn
% Test starts fuly charged so initial soc is assumed as 1
SOC = 1 + cumtrapz(t, I) / 3600 / Q;

%Complete V-t curve
figure('Name','Complete V-t curve','Color','w');
plot(t/3600, V, 'b', 'LineWidth', 1);
grid on; xlabel('Time (h)'); ylabel('Voltage (V)');
title('HPPC test: complete voltage-time curve');

%  FINDING THE PULSES
%  rest_thr tells that if currnet is below 0.5A it is considered as rest
%  "busy" is true while current flows. We pad it with false at both ends so
%  a current block touching the first or last sample is also detected.
%  Each busy block starts at starts(j) and ends at ends(j). Pulses are the
%  short blocks (about 10 s). The rest after a pulse runs from ends(j)+1
%  until the next block starts.

busy   = abs(I) > rest_thr;
b      = [false; busy(:); false];
starts = find(diff(b) ==  1);          % first sample of each current block
ends   = find(diff(b) == -1) - 1;      % last sample of each current block


%  PARAMETERS
%  R0: the voltage jumps by dI*R0 at once when the current starts and again
%       when it stops (ohmic drop, R0 = dV/dI). dI is the REAL current step
%       at each jump (pulse current minus the rest current next to it).
%       We calculate R0 at both jumps and average:
%           R0_start = |V(s)-V(s-1)| / |Ip - I(s-1)|
%           R0_end   = |V(r1)-V(e)|  / |Ip - I(r1)|
%           R0       = (R0_start + R0_end) / 2
%
%  RC: after the current stops the voltage relaxes as
%           V(t) = Vinf + a*exp(-t/RC)
%       We fit this to the rest data over a window of 5RC.
%
%  R1: at the end of the pulse the voltage across the RC branch is
%           dV1 = |a|  =  dI*R1*(1 - exp(-Tp/RC))
%       (the RC branch is NOT fully charged after a short pulse of length
%       Tp, hence the factor in brackets), so
%           R1 = dV1 / ( |dI_start| * (1 - exp(-Tp/RC)) )
%
%  C1: RC = R1*C1

P = struct('type',{},'soc',{},'I',{},'R0',{},'R1',{},'C1',{},'tau',{},'Vinf',{},'a',{},'idx',{});
n = 0;

for j = 1:numel(starts)
    s = starts(j);  e = ends(j);
    dur = t(e) - t(s);
    if dur < pulse_min_dur || dur > pulse_max_dur, continue; end   % not a pulse
    if s < 2 || e >= numel(I), continue; end       % need a sample before and after

    r1 = e + 1;                                    % first rest sample
    if j < numel(starts), r2 = starts(j+1) - 1; else, r2 = numel(I); end

    Ip = median(I(s:e));                           % pulse current [A]
    Tp = t(r1) - t(s);                             % pulse length [s]

    % ---- R0: real current step at each jump, then average ----
    dI_start = Ip - I(s-1);                        % step when the pulse starts
    dI_end   = Ip - I(r1);                         % step when the pulse stops
    R0 = ( abs(V(s) - V(s-1)) / abs(dI_start) ...
         + abs(V(r1) - V(e))  / abs(dI_end) ) / 2;

    % ---- RC from the relaxation fit ----
    [Vinf, a, tau] = fit_relax(t(r1:r2) - t(r1), V(r1:r2), k_rc);

    % ---- R1 and C1 (the RC branch was driven by the step dI_start) ----
    dV1 = abs(a);                                  % RC-branch voltage at pulse end
    R1  = dV1 / (abs(dI_start) * (1 - exp(-Tp/tau)));
    C1  = tau / R1;

    n = n + 1;
    P(n).type = sign(Ip);  P(n).soc = SOC(s);  P(n).I = Ip;
    P(n).R0 = R0;  P(n).R1 = R1;  P(n).C1 = C1;  P(n).tau = tau;
    P(n).Vinf = Vinf;  P(n).a = a;  P(n).idx = [s e r1 r2];
end

% LOOK-UP TABLE HAVING CHARGE AND DISCHARGED PARAMEETRS WITH SOC 
% INTERPOLATED AND EXTRAPOLATED TO COVER WHOLE RANGE
% Between the pulse SOC range we use pchip (smooth). 
% Outside it the values are extrapolated along a straight line.

ptype = [P.type];  psoc = [P.soc];
pR0 = [P.R0];  pR1 = [P.R1];  pC1 = [P.C1];
isD = ptype < 0;   isC = ptype > 0;

SOC_grid = linspace(0, 1, 1001).';                 

R0_charge    = interp_full(psoc(isC), pR0(isC), SOC_grid);   
R1_charge    = interp_full(psoc(isC), pR1(isC), SOC_grid);   
C1_charge    = interp_full(psoc(isC), pC1(isC), SOC_grid);   

R0_discharge = interp_full(psoc(isD), pR0(isD), SOC_grid);   
R1_discharge = interp_full(psoc(isD), pR1(isD), SOC_grid);   
C1_discharge = interp_full(psoc(isD), pC1(isD), SOC_grid);   

LUT = table(SOC_grid, R0_charge, R1_charge, C1_charge, R0_discharge, R1_discharge, C1_discharge,'VariableNames', {'SOC_frac','R0_charge_Ohm','R1_charge_Ohm',    'C1_charge_F','R0_discharge_Ohm', 'R1_discharge_Ohm', 'C1_discharge_F'});

save(fullfile(fileparts(which('setup_paths')),'Lookup','HPPC_params.mat'), 'LUT');

%  FIGURE : R0, R1 and C1 vs SOC (charge and discharge)
x = SOC_grid * 100;
soc_pl = [min(psoc), max(psoc)] * 100;

figure('Name','R0, R1, C1 vs SOC','Color','w');

subplot(3,1,1);
plot(x, R0_charge*1e3, 'r', x, R0_discharge*1e3, 'b', 'LineWidth', 1.5); hold on;
plot(psoc(isC)*100, pR0(isC)*1e3, 'ro', 'MarkerSize', 4);
plot(psoc(isD)*100, pR0(isD)*1e3, 'bo', 'MarkerSize', 4);
yl = ylim;
plot(soc_pl([1 1]), yl, 'k:', 'HandleVisibility', 'off');
plot(soc_pl([2 2]), yl, 'k:', 'HandleVisibility', 'off');
grid on; ylabel('R0 (m\Omega)');
title('HPPC parameters vs SOC  (lines = table, circles = pulses, dotted = pulse SOC range)');
legend('Charge','Discharge','Charge pulses','Discharge pulses','Location','best');

subplot(3,1,2);
plot(x, R1_charge*1e3, 'r', x, R1_discharge*1e3, 'b', 'LineWidth', 1.5); hold on;
plot(psoc(isC)*100, pR1(isC)*1e3, 'ro', 'MarkerSize', 4);
plot(psoc(isD)*100, pR1(isD)*1e3, 'bo', 'MarkerSize', 4);
yl = ylim;
plot(soc_pl([1 1]), yl, 'k:');  plot(soc_pl([2 2]), yl, 'k:');
grid on; ylabel('R1 (m\Omega)');

subplot(3,1,3);
plot(x, C1_charge, 'r', x, C1_discharge, 'b', 'LineWidth', 1.5); hold on;
plot(psoc(isC)*100, pC1(isC), 'ro', 'MarkerSize', 4);
plot(psoc(isD)*100, pC1(isD), 'bo', 'MarkerSize', 4);
yl = ylim;
plot(soc_pl([1 1]), yl, 'k:');  plot(soc_pl([2 2]), yl, 'k:');
grid on; ylabel('C1 (F)'); xlabel('SOC (%)');

% Print the LOOK-UP Table
fprintf('\n   HPPC LOOK-UP TABLE \n');
fprintf('%7s | %10s %10s %10s | %10s %10s %10s\n', 'SOC[%]','R0_chg[Ohm]', 'R1_chg[Ohm]', 'C1_chg[F]','R0_dis[Ohm]', 'R1_dis[Ohm]', 'C1_dis[F]');
for k = 1:100:height(LUT)
    fprintf('%7.1f | %10.5f %10.5f %10.1f | %10.5f %10.5f %10.1f\n', LUT.SOC_frac(k)*100, LUT.R0_charge_Ohm(k), LUT.R1_charge_Ohm(k), LUT.C1_charge_F(k), LUT.R0_discharge_Ohm(k), LUT.R1_discharge_Ohm(k), LUT.C1_discharge_F(k));
end

% FUNCTIONS

% After the pulse is over voltage don't suddnely fall to 0 so that
% behaviour is replicated with help of this function.
% It fits the curve v(t)= Vinf + a*exp(-t/RC)  over a window of k_rc*RC
% Vinf : settled voltage [V]
% a : V(just after the pulse) - Vinf [V]

function [Vinf, a, tau] = fit_relax(tr, vr, k_rc)
    ref = vr(end);
    y   = (vr - ref) * 1e3;                            % mV relative to last sample
    mdl = @(p, x) p(1) + p(2) * exp(-x ./ exp(p(3)));  % p = [Vinf_mV, a_mV, log(RC)]
    sse = @(p, x, yy) sum((mdl(p, x) - yy).^2);
    opt = optimset('Display','off','TolX',1e-10,'TolFun',1e-12,'MaxIter',5000,'MaxFunEvals',5000);

    p = fminsearch(@(q) sse(q, tr, y), [0, y(1), log(30)], opt);   % whole rest

    for it = 1:3                                       % refit on the 5RC window
        mask = tr <= k_rc * exp(p(3));
        if nnz(mask) < 10, mask = true(size(tr)); end
        p = fminsearch(@(q) sse(q, tr(mask), y(mask)), p, opt);
    end

    tau  = exp(p(3));                                  
    Vinf = ref + p(1) / 1e3;                          
    a    = p(2) / 1e3;                                 
end

% INTERPOLATION is done between the pulses and EXTRAPOLATION outside first
% and last pulse over the set.
% Can use a linear interpolation or pchip(smooth) that fits a smooth curve
function yq = interp_full(x, y, xq)
    [x, ix] = sort(x(:));   y = y(:);   y = y(ix);
    yq = interp1(x, y, xq, 'pchip');
    lo = xq < x(1);   hi = xq > x(end);
    yq(lo) = y(1)   + (xq(lo) - x(1))   * (y(2)   - y(1))     / (x(2)   - x(1));
    yq(hi) = y(end) + (xq(hi) - x(end)) * (y(end) - y(end-1)) / (x(end) - x(end-1));
    yq = max(yq, 0.1*min(y));
end