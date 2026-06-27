% hoosier_r20_tire_model.m
%
% fits the pacejka magic formula (bakker-nyborg-pacejka 1989) to fsae ttc round 9 data
% tire: hoosier 43075 16x7.5-10 r20
%
% model predicts lateral force FY as a function of:
%   SA    — slip angle [deg]
%   FZ    — normal force [N]
%   IA    — inclination (camber) angle [deg]
%
% outputs:
%   tireParams struct with fitted pacejka coefficients
%   hoosier_r20_tire_params.mat — saved model file for use in lap sim

clear; clc; close all;

% ------------------------------------------
% section 1: load all hoosier 43075 16x7.5-10 r20 cornering runs


% path to run data folder — update if running from a different working directory
% data now lives in the "Data" subfolder
run_dir = fullfile('Data', 'RunData_Cornering_Matlab_SI_Round9');

all_files = dir(fullfile(run_dir, 'B2356run*.mat'));

fprintf('loading hoosier 43075 16x7.5-10 r20 (8 inch rim) runs...\n');

% preallocate cell arrays to max possible size (one cell per file)
% avoids growing arrays inside the loop — concatenate once after
n_files = length(all_files);
SA_cell = cell(n_files, 1);
FY_cell = cell(n_files, 1);
FZ_cell = cell(n_files, 1);
IA_cell = cell(n_files, 1);
P_cell  = cell(n_files, 1);   % inflation pressure [kPa] — used by optional ~12 psi filter
n_loaded = 0;

for i = 1:n_files
    fpath = fullfile(run_dir, all_files(i).name);

    % load ONLY the tire id first and filter — only hoosier 43075 16x7.5-10 r20.
    % each run file holds ~24 channels of ~60k samples; reading just tireid
    % lets us skip the heavy data entirely for the many non-matching tires.
    % char() handles both cell array and plain char storage formats.
    info = load(fpath, 'tireid');
    tire_str = strtrim(char(info.tireid));
    % only the hoosier 43075 16x7.5-10 r20 on the 8 inch rim.
    % the TTC tested this tire on both 7" and 8" rims; we want 8" only.
    if ~(contains(tire_str, '43075') && contains(tire_str, '16x7.5') ...
            && contains(tire_str, '8 inch rim'))
        continue;
    end

    % matching tire — read the four model channels plus pressure (P, for the
    % optional ~12 psi filter below; loading one extra channel is cheap)
    d = load(fpath, 'SA', 'FY', 'FZ', 'IA', 'P');
    n_loaded = n_loaded + 1;
    SA_cell{n_loaded} = d.SA;
    FY_cell{n_loaded} = d.FY;
    FZ_cell{n_loaded} = d.FZ;
    IA_cell{n_loaded} = d.IA;
    P_cell{n_loaded}  = d.P;
    fprintf('  loaded: %s  (%s)\n', all_files(i).name, tire_str);
end

% fail clearly if the filter matched nothing (wrong folder, or no target tire)
if n_loaded == 0
    error(['no Hoosier 43075 16x7.5-10 R20 (8 inch rim) runs found in "%s". ', ...
           'check that run_dir points at the cornering RunData folder.'], run_dir);
end

% single concatenation after the loop
SA_all = vertcat(SA_cell{1:n_loaded});
FY_all = vertcat(FY_cell{1:n_loaded});
FZ_all = vertcat(FZ_cell{1:n_loaded});
IA_all = vertcat(IA_cell{1:n_loaded});
P_all  = vertcat(P_cell{1:n_loaded});

fprintf('total raw data points: %d\n\n', length(SA_all));

% ----------------------------------------------------------------
% filter out bad data points

% fz is negative in ttc data so take absolute value for model input
FZ_pos = abs(FZ_all);

% discard points where tire is unloaded
fz_min = 100;   % [N]
fz_max = 1300;  % [N] upper physical limit for this tire (300lbs ~= 1300N according to google)

% crop slip angle to valid range
sa_max = 14;    % [deg]

% NOTE: rim width is now fixed to 8" in the load loop above, but matching runs
% still mix inflation pressures (~8-14 psi). that is fine for a single
% pressure-averaged model, but it widens the data cloud. for a tighter, more
% physical fit, additionally filter to your run pressure (the P channel is
% loaded above and the ~12 psi filter below is enabled).
valid = FZ_pos > fz_min & FZ_pos < fz_max & abs(SA_all) < sa_max;

% comment out next line to restrict to 12psi (~82.73 kP)
p_target = 83;  p_tol = 5;                         % [kPa]
valid = valid & abs(P_all - p_target) < p_tol;

SA = SA_all(valid);
FY = FY_all(valid);
FZ = FZ_pos(valid);
IA = IA_all(valid);

% decimate for speed. the TTC logs at >100 Hz, so consecutive samples are
% nearly identical and the dataset is hugely redundant. keeping every Nth
% point makes the nonlinear fit (the slowest step) and the scatter plots run
% several times faster with a negligible effect on the fitted coefficients —
% verified that decim = 5 reproduces the all-data fit to the same rms and
% coefficients. set decim = 1 to fit on every sample.
decim = 3;
SA = SA(1:decim:end);
FY = FY(1:decim:end);
FZ = FZ(1:decim:end);
IA = IA(1:decim:end);

fprintf('data points after filtering: %d  (using every %d-th -> %d points)\n\n', ...
        sum(valid), decim, numel(SA));

% -------------------------------------------------------------------------
% section 3: define pacejka bnp 1989 magic formula structure
%
%   C   = a0                                       (shape factor)
%   D   = a1*Fz^2 + a2*Fz                          (peak force — load sensitivity)
%   BCD = a3*sin(2*atan(Fz/a4))*(1 - a5*|gamma|)   (cornering stiffness; negative in sae convention)
%   E   = a6*Fz + a7                               (curvature factor)
%   SH  = a8*gamma + a9*Fz + a10                   (horizontal shift — plysteer/camber)
%   SV  = a11*Fz*gamma + a12*Fz + a13              (vertical shift — conicity/camber thrust)
%
%   phi = alpha + SH
%   FY  = D*sin(C*atan(B*phi - E*(B*phi - atan(B*phi)))) + SV
%
% note: B = BCD / (C*D)
%       slip angle (alpha) and camber (gamma) in degrees, forces in newtons
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% section 4: fit coefficients (two-stage nonlinear least squares)
%
% This replaces a naive single 14-parameter fit that disagreed with the data.
% Defects that were fixed, and how:
%
%  (A) PEAK-FORCE / LOAD-SENSITIVITY BUG (the important one).
%      Peak lateral force is D = a1*Fz^2 + a2*Fz. With the SAE sign convention
%      D < 0, the normalized peak grip is mu = |D|/Fz = -(a1*Fz + a2). mu only
%      DECREASES with load (real tire load sensitivity) when a1 > 0. An earlier
%      version bounded a1 <= 0, which forced constant mu and made the high-load
%      curves overshoot the data. The old "D concave => a1 <= 0" comment had the
%      sign backwards.
%      Fix: STAGE 1 estimates peak |FY| in load bins directly from the data and
%      least-squares-fits D(Fz) to it, pinning a1, a2 to data-consistent values.
%      A global fit otherwise lets the dense near-zero-slip points wash load
%      sensitivity out, even with a1 freed.
%
%  (B) SHAPE BUG. a0 (C) was sliding to 0.5; below 1 the magic formula cannot
%      form a peak at all. C is now bounded to a physical lateral range
%      [1.2, 1.8]. The curvature offset a7 (E) is allowed positive: the TTC test
%      only sweeps slip to ~13 deg, so the measured curve flattens onto a plateau
%      (E ~ 1) rather than falling off (old bound forced E < 0).
%
%  (C) ROBUSTNESS. STAGE 2 uses a load-normalized residual ((meas-model)/Fz, so
%      every load level counts equally) and a multi-start loop (the fit has many
%      local minima; one start lands in degenerate bound-pinned corners).
%
% NOTE: expect rms ~200 N. That is the raw scatter floor of the TTC data
% (thermal drift, transients, run-to-run spread), not a bad fit. The win here is
% physical coefficients and trustworthy behaviour across load, not a lower rms.
% -------------------------------------------------------------------------

% ---- stage 1: anchor the peak-force load curve D(Fz) to the data ----------
% bin the data by normal force; in each bin estimate peak |FY| as the 90th
% percentile of |FY| in the saturated region (|SA| > 5 deg) near zero camber,
% then fit  -peak = a1*Fz^2 + a2*Fz  (D is negative) by linear least squares.
% (manual percentile via sort -> no Statistics Toolbox dependency.)
fz_anchor_bins = 150:50:1250;          % [N] bin centres
fz_bin_tol     = 55;                   % [N] half-width
fzc = []; Dpk = [];
for fz = fz_anchor_bins
    sel = abs(FZ - fz) < fz_bin_tol & abs(SA) > 5 & abs(IA) < 0.6;
    if sum(sel) > 40
        v = sort(abs(FY(sel)));
        fzc(end+1) = fz;                              %#ok<SAGROW>
        Dpk(end+1) = v(max(1, round(0.90 * numel(v)))); %#ok<SAGROW> peak |FY| proxy
    end
end
fzc = fzc(:); Dpk = Dpk(:);

% solve [Fz^2, Fz] * [a1; a2] = -Dpk  (least squares; D < 0)
ab     = [fzc.^2, fzc] \ (-Dpk);
a1_fix = ab(1);
a2_fix = ab(2);
fprintf('stage 1 — anchored load curve: a1 = %.4e , a2 = %.4f\n', a1_fix, a2_fix);
fprintf('          peak mu  %.2f (at 200 N) -> %.2f (at 1200 N)\n', ...
        abs(a1_fix*200 + a2_fix), abs(a1_fix*1200 + a2_fix));

% ---- stage 2: multi-start, load-normalized fit of the shape parameters ----
% a1, a2 are held at their stage-1 values by setting lb = ub for those two.
%       a0     a1      a2      a3     a4    a5      a6     a7    a8    a9     a10   a11   a12    a13
lb  = [1.2, a1_fix, a2_fix, -5000,   100,   0,  -5e-4,  -1.0,  -2, -1e-2,   -5,  -10,   -2,  -200];
ub  = [1.8, a1_fix, a2_fix,     0,  5000, 0.5,   5e-4,   1.5,   2,  1e-2,    5,   10,    2,   200];

% nominal physically-motivated start (start #1; the rest are random)
%        a0     a1      a2      a3    a4     a5     a6     a7    a8    a9    a10  a11  a12  a13
p_nom = [1.45, a1_fix, a2_fix, -800, 1900,  0.01,  0.0,  0.60, 0.0, 1e-5, 0.0, 0.0, 0.0, 0.0];

% load-normalized residual: (measured - model) / Fz
res_fn = @(p) (FY - mf_lateral(p, SA, FZ, IA)) ./ FZ;

% solver options for lsqnonlin (optimization toolbox required)
opts = optimoptions('lsqnonlin', ...
    'MaxFunctionEvaluations', 60000, ...
    'MaxIterations',          5000,  ...
    'FunctionTolerance',      1e-9,  ...
    'Display',                'off');

n_starts  = 8;
rng(0);                                   % reproducible random starts
best_cost = inf;
p_fit     = p_nom;

fprintf('stage 2 — fitting shape parameters (%d starts)...\n', n_starts);
for s = 1:n_starts
    if s == 1
        p0 = p_nom;                       % deterministic physical start
    else
        p0 = lb + rand(1, numel(lb)) .* (ub - lb);   % random start within bounds
                                                      % (a1,a2 fixed: lb=ub => no spread)
    end

    % guard each start: a pathological random guess (e.g. D ~ 0 -> B blows up)
    % should not abort the whole sweep. the deterministic start #1 guarantees
    % at least one valid fit.
    try
        [p_try, resnorm_try] = lsqnonlin(res_fn, p0, lb, ub, opts);
    catch err
        fprintf('  start %d/%d  skipped (%s)\n', s, n_starts, err.message);
        continue;
    end
    fprintf('  start %d/%d  normalized cost = %.5f\n', s, n_starts, resnorm_try);

    if resnorm_try < best_cost
        best_cost = resnorm_try;
        p_fit     = p_try;
    end
end

% report rms in physical newtons (the optimised cost above is load-normalized)
rms_err = sqrt(mean((FY - mf_lateral(p_fit, SA, FZ, IA)).^2));
fprintf('\nbest normalized cost = %.5f\n', best_cost);
fprintf('fit complete — rms error: %.2f N\n\n', rms_err);

% -------------------------------------------------------------------------
% section 5: store results in a struct for use in lap sim or other tools
% -------------------------------------------------------------------------

tireParams.tire        = 'Hoosier 43075 16x7.5-10 R20';
tireParams.source      = 'FSAE TTC Round 9, Calspan Tire Research Facility';
tireParams.model       = 'Pacejka BNP 1989 Lateral';
tireParams.coeffs      = p_fit;
tireParams.coeff_names = {'a0 (C)', 'a1', 'a2', 'a3', 'a4', 'a5', ...
                           'a6',     'a7', 'a8', 'a9', 'a10', ...
                           'a11',    'a12', 'a13'};
tireParams.rms_error_N = rms_err;

% print coefficient table
fprintf('fitted coefficients:\n');
for k = 1:length(p_fit)
    fprintf('  %-14s = %12.6f\n', tireParams.coeff_names{k}, p_fit(k));
end

% -------------------------------------------------------------------------
% section 6: plot 1 — lateral force vs slip angle at discrete fz levels
%
% solid lines = model, faded scatter = measured data
% primary validation plot: shows how well the model captures measured behavior
% -------------------------------------------------------------------------

fz_bins = [200, 400, 600, 800, 1000, 1200];  % [N] representative fsae corner loads
fz_tol  = 60;                                 % [N] tolerance window around each bin
sa_vec  = linspace(-14, 14, 300);             % slip angle sweep for model curves
cmap    = jet(length(fz_bins));

figure('Name', 'FY vs SA — Hoosier 16x7.5-10 R20', 'Position', [80 80 900 600]);
hold on;

% preallocate to max possible size — trim unused entries after the loop
h_list     = gobjects(length(fz_bins), 1);
label_list = cell(length(fz_bins), 1);
n_plotted  = 0;

for k = 1:length(fz_bins)
    % select data near this fz bin at near-zero camber
    idx = abs(FZ - fz_bins(k)) < fz_tol & abs(IA) < 0.5;
    if sum(idx) < 30; continue; end

    % measured scatter — semi-transparent to reduce visual noise
    scatter(SA(idx), FY(idx), 5, cmap(k,:), 'filled', ...
            'MarkerFaceAlpha', 0.20, 'HandleVisibility', 'off');

    % model curve at this fz level and zero camber
    FY_pred = mf_lateral(p_fit, sa_vec, ...
                          fz_bins(k) * ones(size(sa_vec)), ...
                          zeros(size(sa_vec)));
    h = plot(sa_vec, FY_pred, 'Color', cmap(k,:), 'LineWidth', 2);

    n_plotted = n_plotted + 1;
    h_list(n_plotted)     = h;
    label_list{n_plotted} = sprintf('F_Z = %d N', fz_bins(k));
end

% trim preallocated arrays to the number of bins that had sufficient data
h_list     = h_list(1:n_plotted);
label_list = label_list(1:n_plotted);

xline(0, '--k', 'Alpha', 0.25, 'HandleVisibility', 'off');
yline(0, '--k', 'Alpha', 0.25, 'HandleVisibility', 'off');
xlabel('slip angle  \alpha  [deg]');
ylabel('lateral force  F_Y  [N]');
title('lateral force vs slip angle — IA = 0°  (lines = model, dots = measured)');
legend(h_list, label_list, 'Location', 'best');
grid on; box on;

% -------------------------------------------------------------------------
% section 7: plot 2 — normalized lateral force (mu_y) vs slip angle
%
% mu_y = FY/FZ captures how efficient the tire is at generating lateral grip
% the downward shift of peak mu_y with increasing fz is called load sensitivity —
% a critical input for suspension setup (stiffer arb = more load transfer = less total grip)
% -------------------------------------------------------------------------

figure('Name', 'mu_y vs SA — Load Sensitivity', 'Position', [130 130 900 600]);
hold on;

for k = 1:length(fz_bins)
    FY_pred = mf_lateral(p_fit, sa_vec, ...
                          fz_bins(k) * ones(size(sa_vec)), ...
                          zeros(size(sa_vec)));
    % normalize by normal force to get friction coefficient
    mu_y = FY_pred / fz_bins(k);
    plot(sa_vec, mu_y, 'Color', cmap(k,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('F_Z = %d N', fz_bins(k)));
end

xline(0, '--k', 'Alpha', 0.25, 'HandleVisibility', 'off');
xlabel('slip angle  \alpha  [deg]');
ylabel('\mu_y = F_Y / F_Z  [—]');
title('normalized lateral force vs slip angle — load sensitivity');
legend('Location', 'best');
grid on; box on;

% -------------------------------------------------------------------------
% section 8: plot 3 — camber (ia) effect at a fixed normal force
%
% positive camber (tire leaning inward toward car centerline) shifts and
% increases the lateral force curve — this informs static camber selection
% and camber gain targets in suspension kinematics design
% -------------------------------------------------------------------------

ia_vals  = [0, 2, 4];   % camber angles to compare [deg]
fz_fixed = 600;          % representative mid-corner load for an fsae car [N]

figure('Name', 'Camber Sensitivity — FZ = 600 N', 'Position', [180 180 900 600]);
hold on;
cmap_ia = lines(length(ia_vals));

for k = 1:length(ia_vals)
    FY_pred = mf_lateral(p_fit, sa_vec, ...
                          fz_fixed * ones(size(sa_vec)), ...
                          ia_vals(k) * ones(size(sa_vec)));
    plot(sa_vec, FY_pred, 'Color', cmap_ia(k,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('\\gamma = %d°', ia_vals(k)));
end

xline(0, '--k', 'Alpha', 0.25, 'HandleVisibility', 'off');
xlabel('slip angle  \alpha  [deg]');
ylabel('lateral force  F_Y  [N]');
title(sprintf('camber sensitivity — F_Z = %d N', fz_fixed));
legend('Location', 'best');
grid on; box on;

% -------------------------------------------------------------------------
% section 9: peak grip summary table
%
% peak slip angle at each fz level informs ackermann geometry targets —
% you want the inside/outside wheel slip angles to bracket this value
% -------------------------------------------------------------------------

fprintf('\npeak lateral force summary  (IA = 0 deg):\n');
fprintf('  %-10s  %-14s  %-16s  %-12s\n', ...
        'FZ [N]', '|FY|_peak [N]', 'SA_peak [deg]', 'mu_y_peak');
fprintf('  %s\n', repmat('-', 1, 54));

for k = 1:length(fz_bins)
    FY_curve = mf_lateral(p_fit, sa_vec, ...
                           fz_bins(k) * ones(size(sa_vec)), ...
                           zeros(size(sa_vec)));
    [~, idx_peak] = max(abs(FY_curve));
    peak_fy  = abs(FY_curve(idx_peak));
    peak_sa  = abs(sa_vec(idx_peak));
    peak_mu  = peak_fy / fz_bins(k);
    fprintf('  %-10d  %-14.1f  %-16.2f  %-12.3f\n', ...
            fz_bins(k), peak_fy, peak_sa, peak_mu);
end

% -------------------------------------------------------------------------
% section 10: cornering stiffness vs normal force
%
% cornering stiffness (CS) = dFY/dSA evaluated at SA=0
% from the pacejka model this equals BCD at each FZ level (at zero camber)
% CS is negative in sae convention so we plot |CS| throughout
%
% two subplots:
%   left  — |CS| vs FZ [N/deg]: absolute stiffness, useful for steering torque calcs
%   right — |CS|/FZ vs FZ [1/deg]: normalized stiffness (cornering coefficient)
%           the declining slope here is the linear-range equivalent of load sensitivity —
%           each additional newton of normal force produces less cornering stiffness
%           steeper decline = more penalty for aggressive load transfer
% -------------------------------------------------------------------------

fz_sweep = linspace(100, 1300, 300);   % continuous fz range for model curves [N]

% cornering stiffness from the pacejka bcd term at zero camber
% BCD = a3 * sin(2*atan(Fz/a4)) * (1 - a5*0) = a3 * sin(2*atan(Fz/a4))
CS_model = p_fit(4) .* sin(2 .* atan(fz_sweep ./ p_fit(5)));

% estimate cornering stiffness directly from data using a linear regression
% near SA=0 (where the magic formula is linear): CS = sum(SA*FY)/sum(SA^2)
sa_lin_tol = 2.0;   % [deg] linear range window around SA=0

CS_data  = nan(length(fz_bins), 1);

for k = 1:length(fz_bins)
    idx = abs(FZ - fz_bins(k)) < fz_tol & abs(SA) < sa_lin_tol & abs(IA) < 0.5;
    if sum(idx) < 20; continue; end
    % least-squares slope through origin
    CS_data(k) = (SA(idx)' * FY(idx)) / (SA(idx)' * SA(idx));
end

% use find() to get numeric indices — avoids row/column logical mismatch.
% force both to columns with (:) so the later element-wise divide cs./fz can
% never accidentally broadcast into a matrix, regardless of fz_bins' shape.
idx_valid = find(~isnan(CS_data));
fz_valid  = fz_bins(idx_valid);
fz_valid  = fz_valid(:);               % column vector
cs_valid  = abs(CS_data(idx_valid));
cs_valid  = cs_valid(:);               % column vector

figure('Name', 'Cornering Stiffness vs Normal Force', 'Position', [230 230 1100 520]);

% left subplot: absolute cornering stiffness
subplot(1, 2, 1); hold on;
plot(fz_sweep, abs(CS_model), 'b-', 'LineWidth', 2.5, 'DisplayName', 'model (BCD)');
scatter(fz_valid, cs_valid, 80, 'ro', 'filled', ...
        'DisplayName', 'data (linear fit near SA=0)');
xlabel('normal force  F_Z  [N]');
ylabel('cornering stiffness  |CS|  [N/deg]');
title('cornering stiffness vs normal force');
legend('Location', 'northwest'); grid on; box on;

% right subplot: normalized cornering stiffness
% declining curve shows diminishing returns — each extra N of load buys less stiffness
subplot(1, 2, 2); hold on;
plot(fz_sweep, abs(CS_model) ./ fz_sweep, 'b-', 'LineWidth', 2.5, 'DisplayName', 'model');
scatter(fz_valid, cs_valid ./ fz_valid, 80, 'ro', 'filled', ...
        'DisplayName', 'data estimate');
xlabel('normal force  F_Z  [N]');
ylabel('cornering coefficient  |CS| / F_Z  [1/deg]');
title('normalized cornering stiffness — diminishing returns with load');
legend('Location', 'northeast'); grid on; box on;

sgtitle('hoosier 43075 16x7.5-10 r20  —  cornering stiffness');

% print cornering stiffness table
fprintf('\ncornering stiffness summary  (IA = 0 deg):\n');
fprintf('  %-10s  %-16s  %-20s\n', 'FZ [N]', '|CS| [N/deg]', '|CS|/FZ [1/deg]');
fprintf('  %s\n', repmat('-', 1, 48));
for k = 1:length(fz_bins)
    cs_k = abs(p_fit(4) * sin(2 * atan(fz_bins(k) / p_fit(5))));
    fprintf('  %-10d  %-16.2f  %-20.4f\n', fz_bins(k), cs_k, cs_k / fz_bins(k));
end

% -------------------------------------------------------------------------
% section 10b: plot 4 — 3D lateral force surface  FY(SA, FZ)
%
% the full operating envelope of the fitted 14-coefficient pacejka model:
% lateral force swept across slip angle AND normal force at zero camber.
% the measured points (near IA = 0) are overlaid so the surface can be
% sanity-checked against the raw data it was fit to. change ia_surf below
% to view the surface at a different camber angle.
% -------------------------------------------------------------------------

sa_surf = linspace(-14, 14, 80);        % slip angle grid [deg]
fz_surf = linspace(100, 1300, 80);      % normal force grid [N]
ia_surf = 0;                            % camber for the surface [deg]
[SA_grid, FZ_grid] = meshgrid(sa_surf, fz_surf);

% evaluate the fitted model over the full grid (mf_lateral is element-wise)
FY_surf = mf_lateral(p_fit, SA_grid, FZ_grid, ia_surf * ones(size(SA_grid)));

figure('Name', 'FY surface — Hoosier 16x7.5-10 R20', 'Position', [280 280 950 680]);
surf(SA_grid, FZ_grid, FY_surf, 'EdgeColor', 'none', 'FaceAlpha', 0.90);
hold on;

% overlay measured data near zero camber as a fit sanity check
idx_surf = abs(IA - ia_surf) < 0.5;
scatter3(SA(idx_surf), FZ(idx_surf), FY(idx_surf), 4, ...
         [0.15 0.15 0.15], 'filled', 'MarkerFaceAlpha', 0.15, ...
         'HandleVisibility', 'off');

colormap jet;
cb = colorbar; cb.Label.String = 'lateral force  F_Y  [N]';
xlabel('slip angle  \alpha  [deg]');
ylabel('normal force  F_Z  [N]');
zlabel('lateral force  F_Y  [N]');
title(sprintf(['pacejka lateral force surface — IA = %d°  ', ...
               '(surface = model, dots = measured)'], ia_surf));
view(-135, 25); grid on; box on;

% -------------------------------------------------------------------------
% section 10c: plot 5 — textbook-style lateral force vs slip angle at
%              several camber angles  (the "Fy vs slip angle / camber" figure)
%
% reproduces the classic textbook family of curves: a single representative
% corner load, the slip sweep plotted so FY comes out positive, and one curve
% per inclination angle. for this tire the camber effect is camber THRUST —
% it lifts the low-slip portion of the curve (largest gap near alpha = 0) while
% leaving the peak essentially unchanged, which is exactly the shape shown in
% the reference plot.
%
% the camber sweep uses the TESTED inclination angles (0, 2, 4 deg). the 5/10
% deg values in some textbook figures are outside the TTC test range for this
% tire and would be pure extrapolation, so they are intentionally not used.
% -------------------------------------------------------------------------

fz_ref    = 600;                  % representative fsae corner load [N]
ia_curves = [0, 2, 4];            % tested inclination (camber) angles [deg]
sa_test   = 13;                   % upper slip angle actually swept in the TTC data [deg]
sa_solid  = linspace(0, sa_test, 280);   % within tested range -> solid
sa_dash   = linspace(sa_test, 30, 180);  % beyond tested range -> dashed (extrapolation)
cmap_ia   = lines(length(ia_curves));    % default matlab colors (matches other plots)

figure('Name', 'FY vs SA at multiple camber — textbook style', ...
       'Position', [330 330 900 600]);
hold on;

for k = 1:length(ia_curves)
    % evaluate on the negative-slip branch so lateral force comes out positive,
    % then plot against slip-angle magnitude (matches the "-slip angle"
    % convention of the reference figure)
    FY_solid = mf_lateral(p_fit, -sa_solid, ...
                          fz_ref     * ones(size(sa_solid)), ...
                          ia_curves(k) * ones(size(sa_solid)));
    FY_dash  = mf_lateral(p_fit, -sa_dash, ...
                          fz_ref     * ones(size(sa_dash)), ...
                          ia_curves(k) * ones(size(sa_dash)));

    % solid = fit over measured slip range (this curve appears in the legend)
    plot(sa_solid, FY_solid, 'Color', cmap_ia(k,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('Camber angle = %d°', ia_curves(k)));
    % dashed = extrapolation past the tested range, up to the model peak (~24 deg)
    plot(sa_dash, FY_dash, '--', 'Color', cmap_ia(k,:), 'LineWidth', 2, ...
         'HandleVisibility', 'off');
end

% boundary of the measured data — everything to the right is extrapolated
xline(sa_test, ':k', 'Alpha', 0.4, 'HandleVisibility', 'off');
text(sa_test + 0.3, 150, 'tested range  |  extrapolated', ...
     'Rotation', 90, 'FontSize', 8, 'Color', [0.4 0.4 0.4], ...
     'VerticalAlignment', 'top');

xlabel('slip angle  \alpha  [deg]');
ylabel('lateral force  F_Y  [N]');
title(sprintf('lateral force versus slip angle — F_Z = %d N  (dashed = extrapolated)', fz_ref));
legend('Location', 'southeast');
xlim([0 30]); ylim([0 inf]);
grid on; box on;

% -------------------------------------------------------------------------
% section 11: save tire model parameters
%
% load this file in your lap sim with: load('hoosier_r20_tire_params.mat')
% evaluate the model with: FY = mf_lateral(tireParams.coeffs, SA, FZ, IA)
% -------------------------------------------------------------------------

save('hoosier_r20_tire_params.mat', 'tireParams', 'p_fit');
fprintf('\ntire parameters saved to hoosier_r20_tire_params.mat\n');
fprintf('evaluate model: FY = mf_lateral(tireParams.coeffs, SA, FZ, IA)\n');

% =========================================================================
% local function: pacejka bnp 1989 magic formula — lateral force
%
% call signature: FY = mf_lateral(p, alpha, Fz, gamma)
%   p     — 14-element coefficient vector [a0 ... a13]
%   alpha — slip angle [deg]
%   Fz    — normal force [N], must be positive
%   gamma — inclination (camber) angle [deg]
%   FY    — predicted lateral force [N]
% =========================================================================

function FY = mf_lateral(p, alpha, Fz, gamma)

    % unpack named coefficients from vector
    a0  = p(1);                            % shape factor c
    a1  = p(2);  a2  = p(3);              % peak force d: a1*Fz^2 + a2*Fz
    a3  = p(4);  a4  = p(5);  a5 = p(6); % cornering stiffness bcd
    a6  = p(7);  a7  = p(8);              % curvature factor e
    a8  = p(9);  a9  = p(10); a10 = p(11); % horizontal shift sh
    a11 = p(12); a12 = p(13); a13 = p(14); % vertical shift sv

    % shape factor (typically ~1.3 for lateral)
    C = a0;

    % peak lateral force — quadratic in fz captures nonlinear load sensitivity
    D = a1 .* Fz.^2  +  a2 .* Fz;

    % cornering stiffness — saturates at high load, reduced by camber
    BCD = a3 .* sin(2 .* atan(Fz ./ a4))  .*  (1 - a5 .* abs(gamma));

    % stiffness factor derived from bcd, c, d
    B = BCD ./ (C .* D);

    % curvature factor — controls shape past peak (typically negative or near zero)
    E = a6 .* Fz  +  a7;

    % horizontal shift — camber and load induced offset of the zero-crossing
    SH = a8 .* gamma  +  a9 .* Fz  +  a10;

    % vertical shift — camber thrust and conicity offset
    SV = a11 .* Fz .* gamma  +  a12 .* Fz  +  a13;

    % shifted slip angle
    phi = alpha + SH;

    % magic formula evaluation
    FY = D .* sin( C .* atan( B.*phi - E.*(B.*phi - atan(B.*phi)) ) )  +  SV;

end
