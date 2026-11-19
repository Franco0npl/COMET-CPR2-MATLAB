%% CALIBRATE_COMET_MEASUREMENT_OFFSET_REDUCED_C2
% Dense sweep of a moderate timestamp offset in the reported position signal
% after reducing the artificial servo and physical-feedback degradations.
%
% IMPORTANT: the timestamp offset is not a plant or controller calibration.
% The rebalanced COMET model, controller, actuator dynamics, feedback delay,
% drag and entry ramp are held fixed. Each physical case is simulated once.
% The offset is applied only while recomputing the reported tracking metrics.
%
% A positive measurement lag Delta t means that the position associated with
% timestamp t is x(t - Delta t), while the desired reference remains r(t).
% This tests whether acquisition, estimator or log-synchronisation latency can
% explain part of the experimental position-error gap without changing power.
%
% Fixed rebalanced values:
%   K_cyc = 0.780
%   servo tau = 20 ms
%   physical feedback delay = 2.5 ms
%   horizontal quadratic drag = 0.0525 N/(m/s)^2
%   entry-ramp duration = 1.30 s
%
% Five offsets (40:10:80 ms) x five trajectories x two configurations =
% 50 metric evaluations, but only 10 physical simulations.

clearvars
clc
close all

%% Sweep definition

fixed_k_cyc = 0.780;
fixed_servo_time_constant_s = 0.020;
fixed_feedback_delay_s = 0.0025;
fixed_drag_xy = 0.0525;
fixed_entry_time_s = 1.30;

measurement_lag_candidates_s = [0.040; 0.050; 0.060; 0.070; 0.080];
acceptable_mismatch_percent = 10.0;
save_results = true;

%% Load and verify the selected C2 configuration

script_path = mfilename("fullpath");
script_folder = fileparts(script_path);
addpath(script_folder);

P = comet_parameters();

P.aerodynamics.cyclic_effectiveness_scale = fixed_k_cyc;
P.actuators.servo_time_constant = fixed_servo_time_constant_s;
P.feedback.state_delay = fixed_feedback_delay_s;
P.aerodynamics.quadratic_drag_B(1,1) = fixed_drag_xy;
P.aerodynamics.quadratic_drag_B(2,2) = fixed_drag_xy;

P.disturbance.enabled = false;
P.disturbance.type = "none";
P.trajectory.use_acceleration_limited_entry = true;
P.metrics.evaluate_after_entry = false;

assert(abs(P.aerodynamics.cyclic_effectiveness_scale - fixed_k_cyc) < 1e-12, ...
    "Unexpected K_cyc value.");
assert(abs(P.actuators.servo_time_constant - ...
    fixed_servo_time_constant_s) < 1e-12, ...
    "Unexpected servo time constant.");
assert(abs(P.feedback.state_delay - fixed_feedback_delay_s) < 1e-12, ...
    "Unexpected feedback delay.");
drag_xy = diag(P.aerodynamics.quadratic_drag_B(1:2,1:2));
assert(max(abs(drag_xy - fixed_drag_xy)) < 1e-12, ...
    "Unexpected horizontal quadratic drag.");
assert(all(isfinite(measurement_lag_candidates_s)) && ...
       all(measurement_lag_candidates_s >= 0), ...
    "Measurement-lag candidates must be finite and non-negative.");
assert(numel(unique(measurement_lag_candidates_s)) == ...
       numel(measurement_lag_candidates_s), ...
    "Measurement-lag candidates must be unique.");
assert(P.controller.same_gains_for_all_configurations, ...
    "Single and dual must use identical controller gains.");

tests = P.experimental.tracking;
configurations = ["single", "dual"];

n_lag = numel(measurement_lag_candidates_s);
n_test = numel(tests);
n_configuration = numel(configurations);
n_simulations = n_test*n_configuration;
n_rows = n_lag*n_simulations;

%% Preallocate detailed outputs

Measurement_lag_s = nan(n_rows,1);
Measurement_lag_ms = nan(n_rows,1);
Nominal_shift_at_vmax_m = nan(n_rows,1);
Case = strings(n_rows,1);
Trajectory = strings(n_rows,1);
Configuration = strings(n_rows,1);
Vmax_mps = nan(n_rows,1);
Experimental_completed = false(n_rows,1);

RMSE_apparent_m = nan(n_rows,1);
RMSE_physical_m = nan(n_rows,1);
RMSE_exp_m = nan(n_rows,1);
RMSE_delta_percent = nan(n_rows,1);
Max_error_apparent_m = nan(n_rows,1);
Max_error_physical_m = nan(n_rows,1);
Max_error_exp_m = nan(n_rows,1);
Max_error_delta_percent = nan(n_rows,1);

Average_power_sim_W = nan(n_rows,1);
Average_power_exp_W = nan(n_rows,1);
Average_power_delta_percent = nan(n_rows,1);
Maximum_power_sim_W = nan(n_rows,1);
Maximum_power_exp_W = nan(n_rows,1);
Maximum_power_delta_percent = nan(n_rows,1);

Saturation_percent = nan(n_rows,1);
Peak_cyclic_percent = nan(n_rows,1);
Lost_control = false(n_rows,1);
Tracking_score_percent = nan(n_rows,1);

fprintf("\nCOMET MODERATE MEASUREMENT-OFFSET SWEEP WITH REDUCED FACTORS\n");
fprintf(['Rebalanced model: K_cyc = %.3f | servo tau = %.0f ms | ' ...
    'feedback delay = %.1f ms | d_q,xy = %.4f N/(m/s)^2\n'], ...
    fixed_k_cyc, 1000*fixed_servo_time_constant_s, ...
    1000*fixed_feedback_delay_s, fixed_drag_xy);
fprintf("Fixed entry-ramp duration: %.2f s\n", fixed_entry_time_s);
fprintf("Measurement-lag candidates: %s ms\n", ...
    strjoin(compose("%.0f", 1000*measurement_lag_candidates_s), ", "));
fprintf(['Positive lag compares r(t) with x(t - lag). The controller uses ' ...
    'only the %.1f ms physical feedback delay.\n'], ...
    1000*fixed_feedback_delay_s);
fprintf(['Metrics start at t = 0 | physical simulations: %d | ' ...
    'diagnostic evaluations: %d\n\n'], n_simulations, n_rows);

%% Simulate each physical case once and evaluate every diagnostic offset

row = 0;
simulation_index = 0;
verified_delay_applied_s = NaN;
verified_delay_floor_samples = NaN;
verified_delay_fraction = NaN;

for i_test = 1:n_test
    test = tests(i_test);

    scenario = struct();
    scenario.type = test.trajectory;
    scenario.vmax = test.vmax;
    scenario.amax = test.amax;
    scenario.cycles = P.trajectory.default_cycles;
    scenario.start_from_rest = true;
    scenario.entry_time = fixed_entry_time_s;
    scenario.initial_velocity_I = zeros(3,1);

    for i_configuration = 1:n_configuration
        configuration = configurations(i_configuration);
        target = test.(char(configuration));

        simulation_index = simulation_index + 1;
        fprintf("[%d/%d] %s | %s\n", simulation_index, n_simulations, ...
            test.id, configuration);

        result = simulate_COMET(P, configuration, scenario);
        metrics = result.metrics;

        assert(abs(result.parameters.feedback.state_delay_applied - ...
                   fixed_feedback_delay_s) < 1e-12, ...
            "simulate_COMET did not preserve the fixed feedback delay.");
        assert(abs(result.parameters.aerodynamics.cyclic_effectiveness_scale - ...
                   fixed_k_cyc) < 1e-12, ...
            "simulate_COMET did not preserve the fixed K_cyc.");

        if simulation_index == 1
            verified_delay_applied_s = ...
                result.parameters.feedback.state_delay_applied;
            verified_delay_floor_samples = ...
                result.parameters.feedback.state_delay_floor_samples;
            verified_delay_fraction = ...
                result.parameters.feedback.state_delay_fraction;
        else
            assert(abs(result.parameters.feedback.state_delay_applied - ...
                       verified_delay_applied_s) < 1e-12 && ...
                   result.parameters.feedback.state_delay_floor_samples == ...
                       verified_delay_floor_samples && ...
                   abs(result.parameters.feedback.state_delay_fraction - ...
                       verified_delay_fraction) < 1e-12, ...
                "Feedback-delay implementation changed between cases.");
        end

        has_experiment = target.available && target.completed && ...
            isfinite(target.rmse_position) && ...
            isfinite(target.max_position_error);

        for i_lag = 1:n_lag
            row = row + 1;
            lag_s = measurement_lag_candidates_s(i_lag);

            [rmse_apparent, max_error_apparent] = ...
                apparent_position_metrics(result, lag_s);

            if lag_s == 0
                assert(abs(rmse_apparent - metrics.rmse_position) < 1e-10, ...
                    "Zero-lag RMSE does not reproduce simulate_COMET metrics.");
                assert(abs(max_error_apparent - ...
                    metrics.max_position_error) < 1e-10, ...
                    "Zero-lag maximum error does not reproduce simulate_COMET metrics.");
            end

            Measurement_lag_s(row) = lag_s;
            Measurement_lag_ms(row) = 1000*lag_s;
            Nominal_shift_at_vmax_m(row) = test.vmax*lag_s;
            Case(row) = test.id;
            Trajectory(row) = test.trajectory;
            Configuration(row) = configuration;
            Vmax_mps(row) = test.vmax;
            Experimental_completed(row) = has_experiment;

            RMSE_apparent_m(row) = rmse_apparent;
            RMSE_physical_m(row) = metrics.rmse_position;
            Max_error_apparent_m(row) = max_error_apparent;
            Max_error_physical_m(row) = metrics.max_position_error;

            Average_power_sim_W(row) = metrics.average_power_W;
            Maximum_power_sim_W(row) = metrics.maximum_power_W;
            Saturation_percent(row) = 100*metrics.any_saturation_fraction;
            Peak_cyclic_percent(row) = peak_cyclic_use(result);
            Lost_control(row) = metrics.lost_control;

            if has_experiment
                RMSE_exp_m(row) = target.rmse_position;
                Max_error_exp_m(row) = target.max_position_error;
                Average_power_exp_W(row) = target.average_power;
                Maximum_power_exp_W(row) = target.maximum_power;

                RMSE_delta_percent(row) = percentage_mismatch( ...
                    RMSE_apparent_m(row), RMSE_exp_m(row));
                Max_error_delta_percent(row) = percentage_mismatch( ...
                    Max_error_apparent_m(row), Max_error_exp_m(row));
                Average_power_delta_percent(row) = percentage_mismatch( ...
                    Average_power_sim_W(row), Average_power_exp_W(row));
                Maximum_power_delta_percent(row) = percentage_mismatch( ...
                    Maximum_power_sim_W(row), Maximum_power_exp_W(row));
                Tracking_score_percent(row) = sqrt(mean([ ...
                    RMSE_delta_percent(row)^2, ...
                    Max_error_delta_percent(row)^2]));
            end
        end
    end
end

%% Detailed table

detailed_table = table( ...
    Measurement_lag_s, Measurement_lag_ms, Nominal_shift_at_vmax_m, ...
    Case, Trajectory, Configuration, Vmax_mps, ...
    Experimental_completed, ...
    RMSE_apparent_m, RMSE_physical_m, RMSE_exp_m, ...
    RMSE_delta_percent, ...
    Max_error_apparent_m, Max_error_physical_m, Max_error_exp_m, ...
    Max_error_delta_percent, ...
    Average_power_sim_W, Average_power_exp_W, ...
    Average_power_delta_percent, ...
    Maximum_power_sim_W, Maximum_power_exp_W, ...
    Maximum_power_delta_percent, ...
    Saturation_percent, Peak_cyclic_percent, ...
    Lost_control, Tracking_score_percent);

%% Global diagnostic score

Global_measurement_lag_s = measurement_lag_candidates_s;
Global_measurement_lag_ms = 1000*measurement_lag_candidates_s;
Global_tracking_mismatch_percent = nan(n_lag,1);
Global_power_mismatch_percent = nan(n_lag,1);
Metrics_within_10_percent = nan(n_lag,1);
Power_metrics_within_10_percent = nan(n_lag,1);
Mean_apparent_RMSE_m = nan(n_lag,1);
Mean_apparent_max_error_m = nan(n_lag,1);
Mean_saturation_percent = nan(n_lag,1);
Maximum_saturation_percent = nan(n_lag,1);
Unexpected_loss_count = zeros(n_lag,1);

for i_lag = 1:n_lag
    rows_lag = abs(detailed_table.Measurement_lag_s - ...
        measurement_lag_candidates_s(i_lag)) < 1e-12;
    completed = rows_lag & detailed_table.Experimental_completed;

    mismatch_matrix = [ ...
        detailed_table.RMSE_delta_percent(completed), ...
        detailed_table.Max_error_delta_percent(completed)];
    mismatch_values = mismatch_matrix(isfinite(mismatch_matrix));
    power_mismatch_matrix = [ ...
        detailed_table.Average_power_delta_percent(completed), ...
        detailed_table.Maximum_power_delta_percent(completed)];
    power_mismatch_values = ...
        power_mismatch_matrix(isfinite(power_mismatch_matrix));

    if ~isempty(mismatch_values)
        Global_tracking_mismatch_percent(i_lag) = ...
            sqrt(mean(mismatch_values.^2));
        Metrics_within_10_percent(i_lag) = ...
            100*mean(abs(mismatch_values) <= acceptable_mismatch_percent);
    end

    if ~isempty(power_mismatch_values)
        Global_power_mismatch_percent(i_lag) = ...
            sqrt(mean(power_mismatch_values.^2));
        Power_metrics_within_10_percent(i_lag) = ...
            100*mean(abs(power_mismatch_values) <= ...
            acceptable_mismatch_percent);
    end

    Mean_apparent_RMSE_m(i_lag) = ...
        mean(detailed_table.RMSE_apparent_m(rows_lag));
    Mean_apparent_max_error_m(i_lag) = ...
        mean(detailed_table.Max_error_apparent_m(rows_lag));
    Mean_saturation_percent(i_lag) = ...
        mean(detailed_table.Saturation_percent(rows_lag));
    Maximum_saturation_percent(i_lag) = ...
        max(detailed_table.Saturation_percent(rows_lag));
    Unexpected_loss_count(i_lag) = ...
        sum(detailed_table.Lost_control(rows_lag));
end

global_table = table( ...
    Global_measurement_lag_s, Global_measurement_lag_ms, ...
    Global_tracking_mismatch_percent, Global_power_mismatch_percent, ...
    Metrics_within_10_percent, Power_metrics_within_10_percent, ...
    Mean_apparent_RMSE_m, Mean_apparent_max_error_m, ...
    Mean_saturation_percent, Maximum_saturation_percent, ...
    Unexpected_loss_count);

valid_rows = isfinite(global_table.Global_tracking_mismatch_percent);
assert(any(valid_rows), "No valid experimental comparison is available.");
valid_indices = find(valid_rows);
[~, local_best] = min( ...
    global_table.Global_tracking_mismatch_percent(valid_indices));
selected_global_row = valid_indices(local_best);
selected_lag_s = global_table.Global_measurement_lag_s(selected_global_row);
selected_lag_ms = global_table.Global_measurement_lag_ms(selected_global_row);

%% Print results

fprintf("\nGLOBAL MODERATE MEASUREMENT-OFFSET SCORE\n");
disp(global_table)

fprintf("Selected moderate measurement lag: %.1f ms\n", selected_lag_ms);
fprintf(['Physical feedback-delay check: requested %.1f ms | ' ...
    'applied %.1f ms | floor samples %.0f | fraction %.3f\n'], ...
    1000*fixed_feedback_delay_s, 1000*verified_delay_applied_s, ...
    verified_delay_floor_samples, verified_delay_fraction);
fprintf("Fixed entry-ramp duration: %.2f s\n", fixed_entry_time_s);
fprintf("Global RMSE/max-error mismatch: %.2f %%\n", ...
    global_table.Global_tracking_mismatch_percent(selected_global_row));
fprintf("Experimental tracking metrics within +/-%.0f %%: %.1f %%\n", ...
    acceptable_mismatch_percent, ...
    global_table.Metrics_within_10_percent(selected_global_row));
fprintf("Global average/maximum-power mismatch: %.2f %%\n", ...
    global_table.Global_power_mismatch_percent(selected_global_row));
fprintf(['NOTE: the selected lag changes only reported tracking metrics; ' ...
    'it does not change dynamics, control, power or saturation.\n']);

selected_rows = abs(detailed_table.Measurement_lag_s - ...
    selected_lag_s) < 1e-12;

fprintf("\nCASE-BY-CASE RESULTS AT THE SELECTED DIAGNOSTIC LAG\n");
disp(detailed_table(selected_rows, { ...
    'Case', 'Configuration', 'Experimental_completed', ...
    'Nominal_shift_at_vmax_m', ...
    'RMSE_apparent_m', 'RMSE_physical_m', 'RMSE_exp_m', ...
    'RMSE_delta_percent', ...
    'Max_error_apparent_m', 'Max_error_physical_m', ...
    'Max_error_exp_m', 'Max_error_delta_percent', ...
    'Average_power_sim_W', 'Average_power_exp_W', ...
    'Maximum_power_sim_W', 'Maximum_power_exp_W', ...
    'Saturation_percent', 'Lost_control'}))

%% Results folder

results_folder = fullfile(script_folder, ...
    "measurement_offset_reduced_C2_results");
if save_results && ~isfolder(results_folder)
    mkdir(results_folder)
end

%% Tracking sensitivity figure

figure("Color", "w", ...
    "Name", "COMET moderate measurement-offset sweep", ...
    "Position", [40, 80, 1900, 820]);
tiledlayout(2, n_test, "TileSpacing", "compact", "Padding", "compact");

for i_test = 1:n_test
    nexttile(i_test)
    plot_case_metric(gca, detailed_table, tests(i_test), ...
        "RMSE_apparent_m", "rmse_position", "Position RMSE [m]", ...
        1000*measurement_lag_candidates_s, i_test == 1);
end

for i_test = 1:n_test
    nexttile(n_test + i_test)
    plot_case_metric(gca, detailed_table, tests(i_test), ...
        "Max_error_apparent_m", "max_position_error", ...
        "Maximum position error [m]", ...
        1000*measurement_lag_candidates_s, false);
end

sgtitle(['Moderate timestamp-lag sensitivity with reduced servo and delay; ' ...
    'positive lag compares r(t) with x(t-lag)'])

if save_results
    hide_axes_toolbars(gcf)
    exportgraphics(gcf, fullfile(results_folder, ...
        "COMET_reduced_C2_offset_tracking.png"), "Resolution", 180);
end

%% Global diagnostic figure

figure("Color", "w", "Name", "COMET measurement-offset selection", ...
    "Position", [80, 70, 1450, 820]);
tiledlayout(2,2, "TileSpacing", "compact", "Padding", "compact");

nexttile
plot(global_table.Global_measurement_lag_ms, ...
    global_table.Global_tracking_mismatch_percent, "o-", ...
    "LineWidth", 1.6, "MarkerSize", 6);
hold on
xline(selected_lag_ms, "--", ...
    sprintf("selected %.0f ms", selected_lag_ms));
yline(acceptable_mismatch_percent, "k--", "10% target");
grid on
xlabel("Measurement timestamp lag [ms]")
ylabel("Global tracking mismatch [%]")
title("RMSE and maximum-error score")

nexttile
plot(global_table.Global_measurement_lag_ms, ...
    global_table.Metrics_within_10_percent, "o-", ...
    "LineWidth", 1.6, "MarkerSize", 6);
hold on
xline(selected_lag_ms, "--", ...
    sprintf("selected %.0f ms", selected_lag_ms));
grid on
xlabel("Measurement timestamp lag [ms]")
ylabel("Tracking metrics within +/-10% [%]")
title("Case-level agreement")

nexttile
plot(global_table.Global_measurement_lag_ms, ...
    global_table.Global_power_mismatch_percent, "o-", ...
    "LineWidth", 1.6, "MarkerSize", 6);
hold on
yline(acceptable_mismatch_percent, "k--", "10% target");
grid on
xlabel("Measurement timestamp lag [ms]")
ylabel("Global power mismatch [%]")
title("Power invariance check")

nexttile
plot(global_table.Global_measurement_lag_ms, ...
    global_table.Mean_saturation_percent, "o-", ...
    "LineWidth", 1.5, "DisplayName", "Mean over all cases");
hold on
plot(global_table.Global_measurement_lag_ms, ...
    global_table.Maximum_saturation_percent, "s-", ...
    "LineWidth", 1.5, "DisplayName", "Worst case");
grid on
xlabel("Measurement timestamp lag [ms]")
ylabel("Time in saturation [%]")
title("Physical invariance check")
legend("Location", "best")

if save_results
    hide_axes_toolbars(gcf)
    exportgraphics(gcf, fullfile(results_folder, ...
        "COMET_reduced_C2_offset_selection.png"), "Resolution", 180);
end

%% Selected-lag comparison figure

comparison_rows = selected_rows & detailed_table.Experimental_completed;
selected_comparison = detailed_table(comparison_rows,:);
labels = strrep(selected_comparison.Case, "_", " ") + " " + ...
    upper(extractBetween(selected_comparison.Configuration, 1, 1));

figure("Color", "w", "Name", ...
    "COMET selected measurement-offset comparison", ...
    "Position", [100, 80, 1350, 860]);
tiledlayout(2,2, "TileSpacing", "compact", "Padding", "compact");

nexttile
bar([selected_comparison.RMSE_apparent_m, ...
     selected_comparison.RMSE_exp_m], "grouped")
grid on
ylabel("Position RMSE [m]")
title("Apparent position RMSE")
set(gca, "XTick", 1:height(selected_comparison), ...
    "XTickLabel", labels)
xtickangle(30)
legend("Simulation with timestamp lag", "Experiment", ...
    "Location", "best")

nexttile
bar([selected_comparison.Max_error_apparent_m, ...
     selected_comparison.Max_error_exp_m], "grouped")
grid on
ylabel("Maximum position error [m]")
title("Apparent maximum position error")
set(gca, "XTick", 1:height(selected_comparison), ...
    "XTickLabel", labels)
xtickangle(30)

nexttile
bar([selected_comparison.Average_power_sim_W, ...
     selected_comparison.Average_power_exp_W], "grouped")
grid on
ylabel("Average power [W]")
title("Average power (unchanged)")
set(gca, "XTick", 1:height(selected_comparison), ...
    "XTickLabel", labels)
xtickangle(30)

nexttile
bar([selected_comparison.Maximum_power_sim_W, ...
     selected_comparison.Maximum_power_exp_W], "grouped")
grid on
ylabel("Maximum power [W]")
title("Maximum power (unchanged)")
set(gca, "XTick", 1:height(selected_comparison), ...
    "XTickLabel", labels)
xtickangle(30)

sgtitle(sprintf(['Rebalanced comparison at %.0f ms measurement lag; ' ...
    'physical feedback delay remains %.1f ms'], ...
    selected_lag_ms, 1000*fixed_feedback_delay_s))

if save_results
    hide_axes_toolbars(gcf)
    exportgraphics(gcf, fullfile(results_folder, ...
        "COMET_reduced_C2_offset_selected_comparison.png"), ...
        "Resolution", 180);
end

%% Save numerical outputs

if save_results
    writetable(detailed_table, fullfile(results_folder, ...
        "COMET_reduced_C2_offset_detailed.csv"));
    writetable(global_table, fullfile(results_folder, ...
        "COMET_reduced_C2_offset_global_score.csv"));
    save(fullfile(results_folder, ...
        "COMET_reduced_C2_offset_results.mat"), ...
        "detailed_table", "global_table", ...
        "selected_lag_s", "selected_lag_ms", ...
        "measurement_lag_candidates_s", ...
        "verified_delay_applied_s", ...
        "verified_delay_floor_samples", ...
        "verified_delay_fraction", ...
        "fixed_k_cyc", "fixed_servo_time_constant_s", ...
        "fixed_feedback_delay_s", "fixed_drag_xy", ...
        "fixed_entry_time_s", "acceptable_mismatch_percent");

    fprintf("\nResults saved in:\n%s\n", results_folder);
end


%% Local functions

function [rmse_position, max_position_error] = ...
    apparent_position_metrics(result, measurement_lag_s)

t = result.time(:);
query_time = max(t - measurement_lag_s, t(1));

% Linear interpolation makes this diagnostic continuous and independent of
% the controller sample time. For t < lag, the pre-run vehicle state is
% assumed equal to the initial state (the trajectory starts from rest).
observed_position = interp1(t, result.state.position_I.', ...
    query_time, "linear").';

evaluation = t >= result.metrics.metric_start_time;
if ~any(evaluation)
    evaluation(:) = true;
end

position_error = result.reference.position_I(:,evaluation) - ...
    observed_position(:,evaluation);
position_error_norm = vecnorm(position_error, 2, 1);

rmse_position = sqrt(mean(sum(position_error.^2,1)));
max_position_error = max(position_error_norm);

end


function value = percentage_mismatch(simulated, experimental)

if ~isfinite(experimental) || abs(experimental) < eps
    value = NaN;
else
    value = 100*(simulated - experimental)/abs(experimental);
end

end


function value = peak_cyclic_use(result)

evaluation = result.time >= result.metrics.metric_start_time;
active_rows = 2 + find(logical(result.cyclic_mask(:)));
cyclic = result.actuator.actual(active_rows, evaluation);
limit = max(abs([result.parameters.actuators.cyclic_min, ...
                 result.parameters.actuators.cyclic_max]));

if isempty(cyclic) || limit <= 0
    value = NaN;
else
    value = 100*max(abs(cyclic), [], "all")/limit;
end

end


function plot_case_metric(ax, T, test, simulation_field, ...
    experiment_field, y_label, lag_candidates_ms, show_legend)

single_sim_colour = [0.0000, 0.4470, 0.7410];
single_exp_colour = [0.3010, 0.7450, 0.9330];
dual_sim_colour = [0.8500, 0.3250, 0.0980];
dual_exp_colour = [0.9290, 0.6940, 0.1250];

hold(ax, "on")
x_limits = [min(lag_candidates_ms), max(lag_candidates_ms)];
simulation_values = T.(simulation_field);

draw_experimental_band(ax, x_limits, ...
    test.single, experiment_field, single_exp_colour);
draw_experimental_band(ax, x_limits, ...
    test.dual, experiment_field, dual_exp_colour);

for configuration = ["single", "dual"]
    rows = T.Case == test.id & T.Configuration == configuration;

    if configuration == "single"
        sim_colour = single_sim_colour;
        exp_colour = single_exp_colour;
        marker = "o";
        sim_name = "Timestamp-lag single";
        exp_name = "Experiment single";
        target = test.single;
    else
        sim_colour = dual_sim_colour;
        exp_colour = dual_exp_colour;
        marker = "s";
        sim_name = "Timestamp-lag dual";
        exp_name = "Experiment dual";
        target = test.dual;
    end

    plot(ax, T.Measurement_lag_ms(rows), simulation_values(rows), ...
        "-", "Color", sim_colour, "Marker", marker, ...
        "LineWidth", 1.5, "MarkerSize", 5, ...
        "DisplayName", sim_name);

    if target.available && target.completed && ...
            isfinite(target.(experiment_field))
        yline(ax, target.(experiment_field), "--", ...
            "Color", exp_colour, "LineWidth", 1.6, ...
            "DisplayName", exp_name);
    end
end

if test.single.available && ~test.single.completed
    text(ax, 0.02, 0.96, "Single experiment: loss of control", ...
        "Units", "normalized", "Color", [0.65, 0.1, 0.1], ...
        "VerticalAlignment", "top", "FontSize", 8);
elseif ~test.single.available
    text(ax, 0.02, 0.96, "Single experiment: not reported", ...
        "Units", "normalized", "Color", [0.65, 0.1, 0.1], ...
        "VerticalAlignment", "top", "FontSize", 8);
end

grid(ax, "on")
x_margin = max(2, 0.03*max(diff(x_limits), 1));
xlim(ax, x_limits + [-x_margin, x_margin])
xlabel(ax, "Measurement timestamp lag [ms]")
ylabel(ax, y_label)
title(ax, strrep(test.id, "_", " "))

if show_legend
    legend(ax, "Location", "best", "FontSize", 7);
end

end


function draw_experimental_band(ax, x_limits, target, ...
    experiment_field, colour)

if ~(target.available && target.completed && ...
        isfinite(target.(experiment_field)))
    return
end

experimental = target.(experiment_field);
lower = 0.90*experimental;
upper = 1.10*experimental;

patch(ax, [x_limits(1), x_limits(2), x_limits(2), x_limits(1)], ...
    [lower, lower, upper, upper], colour, ...
    "FaceAlpha", 0.08, "EdgeColor", "none", ...
    "HandleVisibility", "off");

end


function hide_axes_toolbars(fig)

axes_handles = findall(fig, "Type", "axes");
for i = 1:numel(axes_handles)
    try
        axes_handles(i).Toolbar.Visible = "off";
    catch
        % Older MATLAB releases may not expose an axes Toolbar property.
    end
end

end
