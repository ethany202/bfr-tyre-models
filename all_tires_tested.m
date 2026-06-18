% all_tires_tested.m
%
% Compares EVERY tire tested in the FSAE TTC Round 9 cornering dataset.
%
% For each tire the script:
%   1. loads every matching cornering run and filters to clean data,
%   2. fits the Pacejka BNP-1989 lateral magic formula (same model as
%      hoosier_r20_tire_model.m),
%   3. produces a "lateral force vs slip angle at several camber angles"
%      figure (the textbook-style plot) for the BEST rim of each tire model,
%   4. records peak-grip / stiffness metrics for a side-by-side comparison.
%
% Two deliverables are written to the working directory:
%   all_tires_comparison.csv   — one row per tire+rim configuration
%   (figures)                  — one Fy-vs-slip(camber) plot per tire model
%                                plus one overlay comparison figure
%
% Selection of the "best" tire follows the chosen criterion: PEAK LATERAL
% GRIP, mu_y = |Fy|/Fz, evaluated at a representative FSAE corner load.
% Because the same model is tested on two rim widths, each rim is fit
% separately and only the better-gripping rim is plotted / recommended per
% model ("per model, best rim").
%
% -------------------------------------------------------------------------
% NOTE on run time: this fits ~14 tire/rim configurations, each with a
% multi-start nonlinear least-squares fit. Decimation (DECIM) and the number
% of starts (N_STARTS) below trade accuracy for speed — defaults are tuned to
% run in a couple of minutes. Lower DECIM / raise N_STARTS for a tighter fit.
% -------------------------------------------------------------------------

clear; clc; close all;

% ====== user settings ====================================================
run_dir   = fullfile('Data', 'RunData_Cornering_Matlab_SI_Round9');
FZ_MIN    = 100;     % [N]  discard unloaded points
FZ_MAX    = 1300;    % [N]  upper physical limit (~300 lb)
SA_MAX    = 14;      % [deg] crop slip angle
P_TARGET  = 83;      % [kPa] ~12 psi run pressure
P_TOL     = 5;       % [kPa] pressure window
DECIM     = 4;       % keep every Nth sample (speed; data is hugely redundant)
N_STARTS  = 4;       % multi-start count for the shape fit
FZ_REF    = 600;     % [N]  representative FSAE corner load for plots/metrics
FZ_GRIP   = 650;     % [N]  load at which peak-grip ranking is evaluated
IA_CURVES = [0 2 4]; % [deg] tested camber angles for the camber figure
CSV_NAME  = 'all_tires_comparison.csv';
% =========================================================================

% ------------------------------------------------------------------------
% section 1: group every cornering run by its tire id (tire model + rim)
% ------------------------------------------------------------------------
all_files = dir(fullfile(run_dir, 'B2356run*.mat'));
if isempty(all_files)
    error('no B2356run*.mat files found in "%s". check run_dir.', run_dir);
end

ids   = {};               % unique tireid strings
files_by_id = {};         % cell, each a list of file paths for that id

fprintf('scanning %d cornering runs for tire ids...\n', numel(all_files));
for i = 1:numel(all_files)
    fpath = fullfile(run_dir, all_files(i).name);
    info  = load(fpath, 'tireid');
    tid   = strtrim(char(info.tireid));
    j = find(strcmp(ids, tid), 1);
    if isempty(j)
        ids{end+1}        = tid;            %#ok<SAGROW>
        files_by_id{end+1} = {fpath};       %#ok<SAGROW>
    else
        files_by_id{j}{end+1} = fpath;      %#ok<SAGROW>
    end
end
n_cfg = numel(ids);
fprintf('found %d distinct tire/rim configurations:\n', n_cfg);
for k = 1:n_cfg
    fprintf('  [%2d] %-44s (%d runs)\n', k, ids{k}, numel(files_by_id{k}));
end
fprintf('\n');

% ------------------------------------------------------------------------
% section 2: load, filter, fit and characterise each configuration
% ------------------------------------------------------------------------
% preallocate a struct array of results
T = struct('tireid',{},'model',{},'rim_in',{},'n_pts',{}, ...
           'coeffs',{},'rms',{}, ...
           'mu_250',{},'mu_ref',{},'mu_grip',{},'mu_peak',{}, ...
           'Fy_peak_ref',{},'SA_peak_ref',{},'CS_ref',{});

for k = 1:n_cfg
    tid = ids{k};
    [model_name, rim_in] = parse_tireid(tid);
    fprintf('==== [%d/%d] %s ====\n', k, n_cfg, tid);

    % --- load & concatenate the model channels for this config ---
    SA_c = {}; FY_c = {}; FZ_c = {}; IA_c = {}; P_c = {};
    for f = 1:numel(files_by_id{k})
        d = load(files_by_id{k}{f}, 'SA','FY','FZ','IA','P');
        SA_c{end+1}=d.SA; FY_c{end+1}=d.FY; FZ_c{end+1}=d.FZ; %#ok<SAGROW>
        IA_c{end+1}=d.IA; P_c{end+1}=d.P;                     %#ok<SAGROW>
    end
    SA_all=vertcat(SA_c{:}); FY_all=vertcat(FY_c{:});
    FZ_all=vertcat(FZ_c{:}); IA_all=vertcat(IA_c{:}); P_all=vertcat(P_c{:});

    % --- filter ---
    FZ_pos = abs(FZ_all);                       % Fz is negative in TTC data
    valid  = FZ_pos>FZ_MIN & FZ_pos<FZ_MAX & abs(SA_all)<SA_MAX ...
             & abs(P_all-P_TARGET)<P_TOL;
    SA=SA_all(valid); FY=FY_all(valid); FZ=FZ_pos(valid); IA=IA_all(valid);

    % --- decimate ---
    SA=SA(1:DECIM:end); FY=FY(1:DECIM:end);
    FZ=FZ(1:DECIM:end); IA=IA(1:DECIM:end);
    fprintf('  data points after filtering/decimation: %d\n', numel(SA));

    if numel(SA) < 500
        warning('  too few points for "%s" — skipping fit.', tid);
        continue;
    end

    % --- fit Pacejka (two-stage, see fit_pacejka below) ---
    [p_fit, rms_err] = fit_pacejka(SA, FY, FZ, IA, N_STARTS);
    fprintf('  fit rms error: %.1f N\n', rms_err);

    % --- metrics from the fitted model (IA = 0) ---
    sa_vec = linspace(-SA_MAX, SA_MAX, 300);
    mu_at  = @(fz) max(abs(mf_lateral(p_fit, sa_vec, ...
                       fz*ones(size(sa_vec)), zeros(size(sa_vec))))) / fz;
    mu_250  = mu_at(250);
    mu_ref  = mu_at(FZ_REF);
    mu_grip = mu_at(FZ_GRIP);
    % peak mu over a load sweep (the tire's best normalized grip)
    fz_sweep = linspace(FZ_MIN, FZ_MAX, 60);
    mu_sweep = arrayfun(mu_at, fz_sweep);
    mu_peak  = max(mu_sweep);
    % peak force, slip-of-peak, cornering stiffness at the reference load
    FY_ref   = mf_lateral(p_fit, sa_vec, FZ_REF*ones(size(sa_vec)), zeros(size(sa_vec)));
    [Fy_pk, ip] = max(abs(FY_ref));
    SA_pk    = abs(sa_vec(ip));
    CS_ref   = abs(p_fit(4) * sin(2*atan(FZ_REF/p_fit(5))));   % |BCD| at IA=0 [N/deg]

    % --- store ---
    T(end+1) = struct('tireid',tid,'model',model_name,'rim_in',rim_in, ...
        'n_pts',numel(SA),'coeffs',p_fit,'rms',rms_err, ...
        'mu_250',mu_250,'mu_ref',mu_ref,'mu_grip',mu_grip,'mu_peak',mu_peak, ...
        'Fy_peak_ref',Fy_pk,'SA_peak_ref',SA_pk,'CS_ref',CS_ref);  %#ok<SAGROW>
    fprintf('  mu@%dN = %.3f | peak |Fy|@%dN = %.0f N | CS = %.1f N/deg\n\n', ...
        FZ_GRIP, mu_grip, FZ_REF, Fy_pk, CS_ref);
end

if isempty(T)
    error('no tire configurations were successfully fit.');
end

% ------------------------------------------------------------------------
% section 3: pick the best rim per tire model (criterion: peak grip mu)
% ------------------------------------------------------------------------
models = unique({T.model});
best_idx = zeros(numel(models),1);   % index into T of the best rim per model
for m = 1:numel(models)
    sel = find(strcmp({T.model}, models{m}));
    [~, b] = max([T(sel).mu_grip]);
    best_idx(m) = sel(b);
end
% rank the best-rim configs by peak grip (descending)
[~, order] = sort([T(best_idx).mu_grip], 'descend');
best_ranked = best_idx(order);

% ------------------------------------------------------------------------
% section 4: per-model camber figure (best rim) — the reference-style plot
% ------------------------------------------------------------------------
sa_pos  = linspace(0, SA_MAX, 300);
cmap_ia = lines(numel(IA_CURVES));

for r = 1:numel(best_ranked)
    t = T(best_ranked(r));
    figure('Name', sprintf('Fy vs SA — %s', t.model), ...
           'Position', [60+18*r 60+12*r 880 580]);
    hold on;
    for c = 1:numel(IA_CURVES)
        % evaluate on the negative-slip branch so Fy comes out positive,
        % plotted against slip-angle magnitude (matches the reference figure)
        FY_pos = mf_lateral(t.coeffs, -sa_pos, ...
                            FZ_REF*ones(size(sa_pos)), ...
                            IA_CURVES(c)*ones(size(sa_pos)));
        plot(sa_pos, FY_pos, 'Color', cmap_ia(c,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('Camber angle = %d%c', IA_CURVES(c), char(176)));
    end
    xlabel('slip angle  \alpha  [deg]');
    ylabel('lateral force  F_Y  [N]');
    title(sprintf('%s  (%g" rim) — F_Z = %d N', t.model, t.rim_in, FZ_REF));
    legend('Location','east'); xlim([0 SA_MAX]); ylim([0 inf]);
    grid on; box on;
end

% ------------------------------------------------------------------------
% section 5: overlay comparison figure — all best-rim tires, IA = 0
% ------------------------------------------------------------------------
figure('Name','All tires — Fy vs SA comparison (best rim, IA=0)', ...
       'Position',[120 120 980 640]);
hold on;
cmap_all = lines(numel(best_ranked));
for r = 1:numel(best_ranked)
    t = T(best_ranked(r));
    FY_pos = mf_lateral(t.coeffs, -sa_pos, ...
                        FZ_REF*ones(size(sa_pos)), zeros(size(sa_pos)));
    plot(sa_pos, FY_pos, 'Color', cmap_all(r,:), 'LineWidth', 2.2, ...
         'DisplayName', sprintf('%s (%g") — \\mu=%.2f', t.model, t.rim_in, t.mu_grip));
end
xlabel('slip angle  \alpha  [deg]');
ylabel('lateral force  F_Y  [N]');
title(sprintf('lateral force vs slip angle — all tires, best rim, F_Z = %d N, IA = 0%c', ...
      FZ_REF, char(176)));
legend('Location','southeast'); xlim([0 SA_MAX]); ylim([0 inf]);
grid on; box on;

% ------------------------------------------------------------------------
% section 6: write the comparison CSV (all configurations)
% ------------------------------------------------------------------------
% sort the full table by peak grip so the CSV is ranked top-to-bottom
[~, ord_all] = sort([T.mu_grip], 'descend');
Ts = T(ord_all);

is_best = false(numel(Ts),1);
for i = 1:numel(Ts)
    is_best(i) = any(best_ranked == ord_all(i));
end

tbl = table( ...
    {Ts.model}', [Ts.rim_in]', {Ts.tireid}', [Ts.n_pts]', ...
    round([Ts.mu_250]',3), round([Ts.mu_ref]',3), round([Ts.mu_grip]',3), ...
    round([Ts.mu_peak]',3), round([Ts.Fy_peak_ref]',0), round([Ts.SA_peak_ref]',2), ...
    round([Ts.CS_ref]',1), round([Ts.rms]',1), is_best, ...
    'VariableNames', {'tire_model','rim_in','tireid','n_points', ...
       'mu_at_250N','mu_at_600N','mu_at_650N_GRIP','mu_peak', ...
       'Fy_peak_600N_N','SA_peak_deg','CS_600N_N_per_deg','fit_rms_N','best_rim'});

writetable(tbl, CSV_NAME);
fprintf('comparison table written to %s\n\n', CSV_NAME);
disp(tbl);

% ------------------------------------------------------------------------
% section 7: recommendation (criterion: peak lateral grip mu_y)
% ------------------------------------------------------------------------
fprintf('\n========================= RECOMMENDATION =========================\n');
fprintf('Criterion: peak lateral grip mu_y = |Fy|/Fz at ~%d N corner load.\n', FZ_GRIP);
fprintf('(each tire model represented by its higher-gripping rim width)\n\n');
fprintf('  %-4s  %-30s  %-6s  %-8s  %-9s\n','rank','tire model','rim','mu_grip','peak|Fy|');
fprintf('  %s\n', repmat('-',1,66));
for r = 1:numel(best_ranked)
    t = T(best_ranked(r));
    fprintf('  %-4d  %-30s  %g\"     %-8.3f  %-6.0f N\n', ...
        r, t.model, t.rim_in, t.mu_grip, t.Fy_peak_ref);
end
best = T(best_ranked(1));
fprintf('\n>> Best for peak grip: %s on the %g\" rim (mu = %.3f at %d N).\n', ...
        best.model, best.rim_in, best.mu_grip, FZ_GRIP);
fprintf('   Note: the 13\" tires need 13\" wheels; the best 10\"-wheel Hoosier\n');
fprintf('   option is reported in the ranking above for packaging comparison.\n');
fprintf('==================================================================\n');

% =========================================================================
% ============================ local functions ============================
% =========================================================================

function [model_name, rim_in] = parse_tireid(tid)
% split 'Hoosier 43075 16x7.5-10 R20, 7 inch rim' into model + rim width.
    parts = strsplit(tid, ',');
    model_name = strtrim(parts{1});
    rim_in = NaN;
    tok = regexp(tid, '(\d+(\.\d+)?)\s*inch', 'tokens', 'once');
    if ~isempty(tok); rim_in = str2double(tok{1}); end
end

function [p_fit, rms_err] = fit_pacejka(SA, FY, FZ, IA, n_starts)
% two-stage Pacejka BNP-1989 lateral fit (mirrors hoosier_r20_tire_model.m).
% stage 1 anchors the peak-force load curve D(Fz); stage 2 multi-start fits
% the shape parameters with a load-normalized residual.

    % ---- stage 1: anchor D(Fz) = a1*Fz^2 + a2*Fz to binned peak |FY| ----
    fz_bins = 150:50:1250;  tol = 55;
    fzc = []; Dpk = [];
    for fz = fz_bins
        sel = abs(FZ-fz)<tol & abs(SA)>5 & abs(IA)<0.6;
        if sum(sel) > 40
            v = sort(abs(FY(sel)));
            fzc(end+1) = fz;                                 %#ok<AGROW>
            Dpk(end+1) = v(max(1, round(0.90*numel(v))));    %#ok<AGROW>
        end
    end
    if numel(fzc) >= 2
        ab = [fzc(:).^2, fzc(:)] \ (-Dpk(:));
        a1_fix = ab(1); a2_fix = ab(2);
    else
        % fallback if too few bins: rough constant-mu guess
        a1_fix = 0; a2_fix = -2.0;
    end

    % ---- stage 2: bounded multi-start shape fit ----
    %      a0    a1      a2     a3     a4   a5     a6     a7   a8    a9    a10  a11  a12  a13
    lb = [1.2, a1_fix, a2_fix, -5000, 100, 0,  -5e-4, -1.0, -2, -1e-2,  -5, -10,  -2, -200];
    ub = [1.8, a1_fix, a2_fix,     0,5000, 0.5, 5e-4,  1.5,  2,  1e-2,   5,  10,   2,  200];
    p_nom = [1.45, a1_fix, a2_fix, -800, 1900, 0.01, 0.0, 0.60, 0.0, 1e-5, 0.0, 0.0, 0.0, 0.0];

    res_fn = @(p) (FY - mf_lateral(p, SA, FZ, IA)) ./ FZ;
    opts = optimoptions('lsqnonlin', ...
        'MaxFunctionEvaluations', 60000, 'MaxIterations', 5000, ...
        'FunctionTolerance', 1e-9, 'Display', 'off');

    rng(0);  best_cost = inf;  p_fit = p_nom;
    for s = 1:n_starts
        if s == 1
            p0 = p_nom;
        else
            p0 = lb + rand(1,numel(lb)).*(ub-lb);
        end
        try
            [p_try, cost] = lsqnonlin(res_fn, p0, lb, ub, opts);
        catch
            continue;
        end
        if cost < best_cost
            best_cost = cost;  p_fit = p_try;
        end
    end
    rms_err = sqrt(mean((FY - mf_lateral(p_fit, SA, FZ, IA)).^2));
end

function FY = mf_lateral(p, alpha, Fz, gamma)
% Pacejka BNP-1989 magic formula — lateral force.
%   p — 14-element coefficient vector [a0..a13]
%   alpha — slip angle [deg], Fz — normal force [N] (positive), gamma — camber [deg]
    a0=p(1); a1=p(2); a2=p(3); a3=p(4); a4=p(5); a5=p(6); a6=p(7);
    a7=p(8); a8=p(9); a9=p(10); a10=p(11); a11=p(12); a12=p(13); a13=p(14);

    C   = a0;
    D   = a1.*Fz.^2 + a2.*Fz;
    BCD = a3.*sin(2.*atan(Fz./a4)).*(1 - a5.*abs(gamma));
    B   = BCD ./ (C.*D);
    E   = a6.*Fz + a7;
    SH  = a8.*gamma + a9.*Fz + a10;
    SV  = a11.*Fz.*gamma + a12.*Fz + a13;
    phi = alpha + SH;
    FY  = D.*sin(C.*atan(B.*phi - E.*(B.*phi - atan(B.*phi)))) + SV;
end
