%% CALIBRATE_COMET_DELAY_FIXED_RAMP
% Targeted sweep of the common feedback-state delay with a fixed entry ramp.
%
% The entry ramp is fixed at 1.30 s, the value that reproduced the reported
% average and maximum power without exceeding the physical guardrails. RMSE,
% maximum position error, power and saturation are evaluated from t = 0.
%
% Only the delay applied to the state seen by the controller changes. The
% plant state is never delayed. Fractional delays are reconstructed by the
% interpolation already implemented in simulate_COMET.
%
% Fixed C2 values:
%   K_cyc = 0.780, servo tau = 30 ms,
%   horizontal quadratic drag = 0.0525 N/(m/s)^2,
%   entry-ramp duration = 1.30 s.
%
% Five delays x five trajectories x two configurations = 50 runs.

clearvars
clc
close all

%% Sweep definition

fixed_entry_time_s = 1.30;
delay_candidates_s = [0.0050; 0.0075; 0.0100; 0.0125; 0.0150];
acceptable_mismatch_percent = 10.0;
power_guardrail_percent = 15.0;
saturation_guardrail_percent = 5.0;
save_results = true;

%% Load and verify the selected C2 parameter configuration

script_path = mfilename("fullpath");
script_folder = fileparts(script_path);
addpath(script_folder);

P_base = comet_parameters();

% Apply C2 explicitly so this calibration remains reproducible even when an
% older comet_parameters.m is still present elsewhere on the MATLAB path.
P_base.aerodynamics.cyclic_effectiveness_scale = 0.780;
P_base.actuators.servo_time_constant = 0.030;
P_base.feedback.state_delay = delay_candidates_s(1);
P_base.aerodynamics.quadratic_drag_B(1,1) = 0.0525;
P_base.aerodynamics.quadratic_drag_B(2,2) = 0.0525;

P_base.disturbance.enabled = false;
P_base.disturbance.type = "none";
P_base.trajectory.use_acceleration_limited_entry = true;
P_base.metrics.evaluate_after_entry = false;

assert(abs(P_base.aerodynamics.cyclic_effectiveness_scale - 0.78) < 1e-12, ...
    "This sweep expects K_cyc = 0.78.");
assert(abs(P_base.actuators.servo_time_constant - 0.030) < 1e-12, ...
    "This sweep expects the C2 servo_time_constant = 0.030 s.");
assert(all(isfinite(delay_candidates_s)) && ...
       all(delay_candidates_s >= 0), ...
    "All feedback-delay candidates must be finite and nonnegative.");
assert(numel(unique(delay_candidates_s)) == numel(delay_candidates_s), ...
    "Feedback-delay candidates must be unique.");
drag_xy = diag(P_base.aerodynamics.quadratic_drag_B(1:2,1:2));
assert(max(abs(drag_xy - 0.0525)) < 1e-12, ...
    "This sweep expects the C2 horizontal quadratic drag = 0.0525 N/(m/s)^2.");
assert(P_base.controller.same_gains_for_all_configurations, ...
    "Single and dual must use identical controller gains.");

tests = P_base.experimental.tracking;
configurations = ["single", "dual"];

n_delay = numel(delay_candidates_s);
n_test = numel(tests);
n_configuration = numel(configurations);
n_rows = n_delay*n_test*n_configuration;

%% Preallocate detailed outputs

Delay_requested_s = nan(n_rows,1);
Delay_requested_ms = nan(n_rows,1);
Delay_applied_s = nan(n_rows,1);
Delay_samples = nan(n_rows,1);
Case = strings(n_rows,1);
Trajectory = strings(n_rows,1);
Configuration = strings(n_rows,1);
Vmax_mps = nan(n_rows,1);
Experimental_completed = false(n_rows,1);

RMSE_sim_m = nan(n_rows,1);
RMSE_exp_m = nan(n_rows,1);
RMSE_delta_percent = nan(n_rows,1);
Max_error_sim_m = nan(n_rows,1);
Max_error_exp_m = nan(n_rows,1);
Max_error_delta_percent = nan(n_rows,1);
Velocity_RMSE_sim_mps = nan(n_rows,1);
Average_power_sim_W = nan(n_rows,1);
Average_power_exp_W = nan(n_rows,1);
Average_power_delta_percent = nan(n_rows,1);
Maximum_power_sim_W = nan(n_rows,1);
Maximum_power_exp_W = nan(n_rows,1);
Maximum_power_delta_percent = nan(n_rows,1);

Roll_peak_sim = nan(n_rows,1);
Roll_peak_exp = nan(n_rows,1);
Pitch_peak_sim = nan(n_rows,1);
Pitch_peak_exp = nan(n_rows,1);

Saturation_percent = nan(n_rows,1);
Cyclic_saturation_percent = nan(n_rows,1);
Rate_saturation_percent = nan(n_rows,1);
Peak_cyclic_percent = nan(n_rows,1);
Lost_control = false(n_rows,1);
Tracking_score_percent = nan(n_rows,1);

fprintf("\nCOMET C2 FEEDBACK-DELAY SWEEP WITH FIXED ENTRY RAMP\n");
fprintf("K_cyc = %.3f | servo tau = %.0f ms | d_q,xy = %.4f N/(m/s)^2\n", ...
    P_base.aerodynamics.cyclic_effectiveness_scale, ...
    1000*P_base.actuators.servo_time_constant, ...
    drag_xy(1));
fprintf("Fixed entry-ramp duration: %.2f s\n", fixed_entry_time_s);
fprintf("Feedback delays: %s ms\n", ...
    strjoin(compose("%.1f", 1000*delay_candidates_s), ", "));
fprintf("Metrics start at t = 0 | total simulations: %d\n\n", n_rows);

%% Run every delay, test and configuration

row = 0;
for i_delay = 1:n_delay
    fprintf("[%d/%d] feedback delay = %.1f ms\n", ...
        i_delay, n_delay, 1000*delay_candidates_s(i_delay));

    P = P_base;
    P.feedback.state_delay = delay_candidates_s(i_delay);

    for i_test = 1:n_test
        test = tests(i_test);

        scenario = struct();
        scenario.type = test.trajectory;
        scenario.vmax = test.vmax;
        scenario.amax = test.amax;
        scenario.cycles = P_base.trajectory.default_cycles;
        scenario.start_from_rest = true;
        scenario.entry_time = fixed_entry_time_s;
        scenario.initial_velocity_I = zeros(3,1);

        for i_configuration = 1:n_configuration
            row = row + 1;
            configuration = configurations(i_configuration);
            target = test.(char(configuration));

            result = simulate_COMET(P, configuration, scenario);
            metrics = result.metrics;

            Delay_requested_s(row) = delay_candidates_s(i_delay);
            Delay_requested_ms(row) = 1000*delay_candidates_s(i_delay);
            Delay_applied_s(row) = ...
                result.parameters.feedback.state_delay_applied;
            Delay_samples(row) = ...
                result.parameters.feedback.state_delay_samples;
            Case(row) = test.id;
            Trajectory(row) = test.trajectory;
            Configuration(row) = configuration;
            Vmax_mps(row) = test.vmax;

            RMSE_sim_m(row) = metrics.rmse_position;
            Max_error_sim_m(row) = metrics.max_position_error;
            Velocity_RMSE_sim_mps(row) = metrics.rmse_velocity;
            Average_power_sim_W(row) = metrics.average_power_W;
            Maximum_power_sim_W(row) = metrics.maximum_power_W;
            Roll_peak_sim(row) = max(abs(metrics.roll_command_range));
            Pitch_peak_sim(row) = max(abs(metrics.pitch_command_range));

            Saturation_percent(row) = 100*metrics.any_saturation_fraction;
            Cyclic_saturation_percent(row) = ...
                100*metrics.cyclic_position_saturation_fraction;
            Rate_saturation_percent(row) = ...
                100*metrics.actuator_rate_saturation_fraction;
            Peak_cyclic_percent(row) = peak_cyclic_use(result);
            Lost_control(row) = metrics.lost_control;

            has_experiment = target.available && target.completed && ...
                isfinite(target.rmse_position) && ...
                isfinite(target.max_position_error);
            Experimental_completed(row) = has_experiment;

            if has_experiment
                RMSE_exp_m(row) = target.rmse_position;
                Max_error_exp_m(row) = target.max_position_error;
                Average_power_exp_W(row) = target.average_power;
                Maximum_power_exp_W(row) = target.maximum_power;
                Roll_peak_exp(row) = max(abs(target.roll_command_range));
                Pitch_peak_exp(row) = max(abs(target.pitch_command_range));

                RMSE_delta_percent(row) = percentage_mismatch( ...
                    RMSE_sim_m(row), RMSE_exp_m(row));
                Max_error_delta_percent(row) = percentage_mismatch( ...
                    Max_error_sim_m(row), Max_error_exp_m(row));
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
    Delay_requested_s, Delay_requested_ms, Delay_applied_s, ...
    Delay_samples, Case, Trajectory, Configuration, Vmax_mps, ...
    Experimental_completed, ...
    RMSE_sim_m, RMSE_exp_m, RMSE_delta_percent, ...
    Max_error_sim_m, Max_error_exp_m, Max_error_delta_percent, ...
    Velocity_RMSE_sim_mps, ...
    Average_power_sim_W, Average_power_exp_W, ...
    Average_power_delta_percent, ...
    Maximum_power_sim_W, Maximum_power_exp_W, ...
    Maximum_power_delta_percent, ...
    Roll_peak_sim, Roll_peak_exp, ...
    Pitch_peak_sim, Pitch_peak_exp, ...
    Saturation_percent, Cyclic_saturation_percent, ...
    Rate_saturation_percent, Peak_cyclic_percent, ...
    Lost_control, Tracking_score_percent);

%% Global score and physical guardrails

Global_delay_s = delay_candidates_s;
Global_delay_ms = 1000*delay_candidates_s;
Global_tracking_mismatch_percent = nan(n_delay,1);
Global_power_mismatch_percent = nan(n_delay,1);
Metrics_within_10_percent = nan(n_delay,1);
Power_metrics_within_10_percent = nan(n_delay,1);
Mean_saturation_percent = nan(n_delay,1);
Maximum_saturation_percent = nan(n_delay,1);
Completed_saturation_over_guardrail_count = zeros(n_delay,1);
Unexpected_loss_count = zeros(n_delay,1);

for i_delay = 1:n_delay
    completed = abs(detailed_table.Delay_requested_s - ...
        delay_candidates_s(i_delay)) < 1e-12 & ...
        detailed_table.Experimental_completed;
    all_rows = abs(detailed_table.Delay_requested_s - ...
        delay_candidates_s(i_delay)) < 1e-12;

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
        Global_tracking_mismatch_percent(i_delay) = ...
            sqrt(mean(mismatch_values.^2));
        Metrics_within_10_percent(i_delay) = ...
            100*mean(abs(mismatch_values) <= acceptable_mismatch_percent);
    end

    if ~isempty(power_mismatch_values)
        Global_power_mismatch_percent(i_delay) = ...
            sqrt(mean(power_mismatch_values.^2));
        Power_metrics_within_10_percent(i_delay) = ...
            100*mean(abs(power_mismatch_values) <= ...
            acceptable_mismatch_percent);
    end

    Mean_saturation_percent(i_delay) = ...
        mean(detailed_table.Saturation_percent(all_rows));
    Maximum_saturation_percent(i_delay) = ...
        max(detailed_table.Saturation_percent(all_rows));
    Completed_saturation_over_guardrail_count(i_delay) = sum( ...
        detailed_table.Saturation_percent(completed) > ...
        saturation_guardrail_percent);
    Unexpected_loss_count(i_delay) = sum( ...
        detailed_table.Lost_control(completed));
end

global_table = table( ...
    Global_delay_s, Global_delay_ms, Global_tracking_mismatch_percent, ...
    Global_power_mismatch_percent, Metrics_within_10_percent, ...
    Power_metrics_within_10_percent, Mean_saturation_percent, ...
    Maximum_saturation_percent, ...
    Completed_saturation_over_guardrail_count, ...
    Unexpected_loss_count);

physical_eligible = global_table.Unexpected_loss_count == 0 & ...
    global_table.Completed_saturation_over_guardrail_count == 0 & ...
    isfinite(global_table.Global_tracking_mismatch_percent);

eligible = physical_eligible & ...
    global_table.Global_power_mismatch_percent <= power_guardrail_percent;

if ~any(eligible)
    warning(["No delay satisfies the physical and power guardrails. " ...
        "Selection is restricted to physically eligible delays."]);
    eligible = physical_eligible;
end
if ~any(eligible)
    warning(["Every delay violates at least one physical guardrail. " ...
        "Selection is restricted to delays with no unexpected loss of control."]);
    eligible = global_table.Unexpected_loss_count == 0 & ...
        isfinite(global_table.Global_tracking_mismatch_percent);
end
if ~any(eligible)
    warning("Every delay loses control in at least one completed experiment.");
    eligible = isfinite(global_table.Global_tracking_mismatch_percent);
end

eligible_rows = find(eligible);
[~, local_best] = min( ...
    global_table.Global_tracking_mismatch_percent(eligible_rows));
selected_global_row = eligible_rows(local_best);
selected_delay_s = global_table.Global_delay_s(selected_global_row);
selected_delay_ms = global_table.Global_delay_ms(selected_global_row);

%% Print concise results

fprintf("\nGLOBAL FEEDBACK-DELAY SCORE\n");
disp(global_table)

fprintf("Selected feedback delay: %.1f ms\n", selected_delay_ms);
fprintf("Fixed entry-ramp duration: %.2f s\n", fixed_entry_time_s);
fprintf("Global RMSE/max-error mismatch: %.2f %%\n", ...
    global_table.Global_tracking_mismatch_percent(selected_global_row));
fprintf("Experimental tracking metrics within +/-%.0f %%: %.1f %%\n", ...
    acceptable_mismatch_percent, ...
    global_table.Metrics_within_10_percent(selected_global_row));
fprintf("Global average/maximum-power mismatch: %.2f %%\n", ...
    global_table.Global_power_mismatch_percent(selected_global_row));
fprintf("Experimental power metrics within +/-%.0f %%: %.1f %%\n", ...
    acceptable_mismatch_percent, ...
    global_table.Power_metrics_within_10_percent(selected_global_row));
fprintf("Completed cases above %.1f %% saturation: %d\n", ...
    saturation_guardrail_percent, ...
    global_table.Completed_saturation_over_guardrail_count(selected_global_row));

selected_rows = abs(detailed_table.Delay_requested_s - ...
    selected_delay_s) < 1e-12;

fprintf("\nCASE-BY-CASE RESULTS AT THE SELECTED FEEDBACK DELAY\n");
disp(detailed_table(selected_rows, { ...
    'Case', 'Configuration', 'Experimental_completed', ...
    'RMSE_sim_m', 'RMSE_exp_m', 'RMSE_delta_percent', ...
    'Max_error_sim_m', 'Max_error_exp_m', ...
    'Max_error_delta_percent', ...
    'Average_power_sim_W', 'Average_power_exp_W', ...
    'Average_power_delta_percent', ...
    'Maximum_power_sim_W', 'Maximum_power_exp_W', ...
    'Maximum_power_delta_percent', ...
    'Saturation_percent', 'Peak_cyclic_percent', ...
    'Lost_control'}))

%% Combined tracking figure: one window, both configurations

results_folder = fullfile(script_folder, ...
    "feedback_delay_C2_fixed_ramp_sweep_results");
if save_results && ~isfolder(results_folder)
    mkdir(results_folder)
end

figure("Color", "w", "Name", "COMET feedback-delay tracking sweep", ...
    "Position", [40, 80, 1900, 820]);
tiledlayout(2, n_test, "TileSpacing", "compact", "Padding", "compact");

for i_test = 1:n_test
    nexttile(i_test)
    plot_case_metric(gca, detailed_table, tests(i_test), ...
        "RMSE_sim_m", "rmse_position", "Position RMSE [m]", ...
        1000*delay_candidates_s, saturation_guardrail_percent, ...
        i_test == 1);
end

for i_test = 1:n_test
    nexttile(n_test + i_test)
    plot_case_metric(gca, detailed_table, tests(i_test), ...
        "Max_error_sim_m", "max_position_error", ...
        "Maximum position error [m]", 1000*delay_candidates_s, ...
        saturation_guardrail_percent, false);
end

sgtitle(["C2 feedback-delay sensitivity at 1.30 s entry ramp; " ...
    "metrics evaluated from t = 0; " ...
    "solid = simulation, dashed = experiment"])

if save_results
    hide_axes_toolbars(gcf)
    exportgraphics(gcf, fullfile(results_folder, ...
        "COMET_feedback_delay_tracking.png"), "Resolution", 180);
end

%% Combined power figure: average and maximum power in one window

figure("Color", "w", "Name", "COMET C2 feedback-delay power sweep", ...
    "Position", [40, 80, 1900, 820]);
tiledlayout(2, n_test, "TileSpacing", "compact", "Padding", "compact");

for i_test = 1:n_test
    nexttile(i_test)
    plot_case_metric(gca, detailed_table, tests(i_test), ...
        "Average_power_sim_W", "average_power", "Average power [W]", ...
        1000*delay_candidates_s, saturation_guardrail_percent, ...
        i_test == 1);
end

for i_test = 1:n_test
    nexttile(n_test + i_test)
    plot_case_metric(gca, detailed_table, tests(i_test), ...
        "Maximum_power_sim_W", "maximum_power", "Maximum power [W]", ...
        1000*delay_candidates_s, saturation_guardrail_percent, false);
end

sgtitle(["C2 feedback-delay power sensitivity at 1.30 s entry ramp; " ...
    "metrics evaluated from t = 0; " ...
    "solid = simulation, dashed = experiment"])

if save_results
    hide_axes_toolbars(gcf)
    exportgraphics(gcf, fullfile(results_folder, ...
        "COMET_feedback_delay_power.png"), "Resolution", 180);
end

%% Global score and guardrails figure

figure("Color", "w", "Name", "COMET feedback-delay selection", ...
    "Position", [80, 70, 1450, 820]);
tiledlayout(2,2, "TileSpacing", "compact", "Padding", "compact");

nexttile
plot(global_table.Global_delay_ms, ...
    global_table.Global_tracking_mismatch_percent, "o-", ...
    "LineWidth", 1.6, "MarkerSize", 6);
hold on
xline(selected_delay_ms, "--", ...
    sprintf("selected %.1f ms", selected_delay_ms));
yline(acceptable_mismatch_percent, "k--", "10% target");
grid on
xlabel("Feedback delay [ms]")
ylabel("Global tracking mismatch [%]")
title("RMSE and maximum-error score")

nexttile
plot(global_table.Global_delay_ms, ...
    global_table.Global_power_mismatch_percent, "o-", ...
    "LineWidth", 1.6, "MarkerSize", 6);
hold on
xline(selected_delay_ms, "--", ...
    sprintf("tracking-selected %.1f ms", selected_delay_ms));
yline(acceptable_mismatch_percent, "k--", "10% target");
yline(power_guardrail_percent, "r:", ...
    sprintf("%.0f%% selection guardrail", power_guardrail_percent));
grid on
xlabel("Feedback delay [ms]")
ylabel("Global power mismatch [%]")
title("Average- and maximum-power score")

nexttile
plot(global_table.Global_delay_ms, ...
    global_table.Mean_saturation_percent, "o-", ...
    "LineWidth", 1.5, "DisplayName", "Mean over all cases");
hold on
plot(global_table.Global_delay_ms, ...
    global_table.Maximum_saturation_percent, "s-", ...
    "LineWidth", 1.5, "DisplayName", "Worst case");
yline(saturation_guardrail_percent, "k--", ...
    sprintf("%.0f%% guardrail", saturation_guardrail_percent));
grid on
xlabel("Feedback delay [ms]")
ylabel("Time in saturation [%]")
title("Saturation exposure")
legend("Location", "best")

nexttile
bar(global_table.Global_delay_ms, [ ...
    global_table.Completed_saturation_over_guardrail_count, ...
    global_table.Unexpected_loss_count], "grouped");
grid on
xlabel("Feedback delay [ms]")
ylabel("Number of completed experimental cases")
title("Physical guardrail violations")
legend("Saturation > guardrail", "Unexpected loss", ...
    "Location", "best")

if save_results
    hide_axes_toolbars(gcf)
    exportgraphics(gcf, fullfile(results_folder, ...
        "COMET_feedback_delay_selection.png"), "Resolution", 180);
end

%% Selected-delay simulation-versus-experiment figure

comparison_rows = selected_rows & detailed_table.Experimental_completed;
selected_comparison = detailed_table(comparison_rows,:);
labels = strrep(selected_comparison.Case, "_", " ") + " " + ...
    upper(extractBetween(selected_comparison.Configuration, 1, 1));

figure("Color", "w", "Name", "COMET selected feedback-delay comparison", ...
    "Position", [100, 80, 1350, 860]);
tiledlayout(2,2, "TileSpacing", "compact", "Padding", "compact");

nexttile
bar([selected_comparison.RMSE_sim_m, ...
     selected_comparison.RMSE_exp_m], "grouped")
grid on
ylabel("Position RMSE [m]")
title("Position RMSE")
set(gca, "XTick", 1:height(selected_comparison), ...
    "XTickLabel", labels)
xtickangle(30)
legend("Simulation", "Experiment", "Location", "best")

nexttile
bar([selected_comparison.Max_error_sim_m, ...
     selected_comparison.Max_error_exp_m], "grouped")
grid on
ylabel("Maximum position error [m]")
title("Maximum position error")
set(gca, "XTick", 1:height(selected_comparison), ...
    "XTickLabel", labels)
xtickangle(30)

nexttile
bar([selected_comparison.Average_power_sim_W, ...
     selected_comparison.Average_power_exp_W], "grouped")
grid on
ylabel("Average power [W]")
title("Average power")
set(gca, "XTick", 1:height(selected_comparison), ...
    "XTickLabel", labels)
xtickangle(30)

nexttile
bar([selected_comparison.Maximum_power_sim_W, ...
     selected_comparison.Maximum_power_exp_W], "grouped")
grid on
ylabel("Maximum power [W]")
title("Maximum power")
set(gca, "XTick", 1:height(selected_comparison), ...
    "XTickLabel", labels)
xtickangle(30)

sgtitle(sprintf( ...
    "Simulation versus experiment at %.1f ms delay and %.2f s entry ramp", ...
    selected_delay_ms, fixed_entry_time_s))

if save_results
    hide_axes_toolbars(gcf)
    exportgraphics(gcf, fullfile(results_folder, ...
        "COMET_feedback_delay_selected_comparison.png"), ...
        "Resolution", 180);
end

%% Save numerical outputs

if save_results
    writetable(detailed_table, fullfile(results_folder, ...
        "COMET_feedback_delay_detailed.csv"));
    writetable(global_table, fullfile(results_folder, ...
        "COMET_feedback_delay_global_score.csv"));
    save(fullfile(results_folder, ...
        "COMET_feedback_delay_sweep_results.mat"), ...
        "detailed_table", "global_table", ...
        "selected_delay_s", "selected_delay_ms", ...
        "delay_candidates_s", "fixed_entry_time_s", ...
        "acceptable_mismatch_percent", ...
        "power_guardrail_percent", ...
        "saturation_guardrail_percent");

    fprintf("\nResults saved in:\n%s\n", results_folder);
end


%% Local functions

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
    experiment_field, y_label, delay_candidates_ms, saturation_guardrail, ...
    show_legend)

single_sim_colour = [0.0000, 0.4470, 0.7410];
single_exp_colour = [0.3010, 0.7450, 0.9330];
dual_sim_colour = [0.8500, 0.3250, 0.0980];
dual_exp_colour = [0.9290, 0.6940, 0.1250];

hold(ax, "on")
x_limits = [min(delay_candidates_ms), max(delay_candidates_ms)];
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
        sim_name = "Simulation single";
        exp_name = "Experiment single";
        target = test.single;
    else
        sim_colour = dual_sim_colour;
        exp_colour = dual_exp_colour;
        marker = "s";
        sim_name = "Simulation dual";
        exp_name = "Experiment dual";
        target = test.dual;
    end

    plot(ax, T.Delay_requested_ms(rows), simulation_values(rows), ...
        "-", "Color", sim_colour, "Marker", marker, ...
        "LineWidth", 1.5, "MarkerSize", 5, ...
        "DisplayName", sim_name);

    if target.available && target.completed && ...
            isfinite(target.(experiment_field))
        yline(ax, target.(experiment_field), "--", ...
            "Color", exp_colour, "LineWidth", 1.6, ...
            "DisplayName", exp_name);
    end

    saturated = rows & T.Saturation_percent > saturation_guardrail;
    if any(saturated)
        plot(ax, T.Delay_requested_ms(saturated), ...
            simulation_values(saturated), "kx", ...
            "LineWidth", 1.3, "MarkerSize", 7, ...
            "HandleVisibility", "off");
    end

    lost = rows & T.Lost_control;
    if any(lost)
        plot(ax, T.Delay_requested_ms(lost), ...
            simulation_values(lost), "r*", ...
            "LineWidth", 1.1, "MarkerSize", 8, ...
            "HandleVisibility", "off");
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
x_margin = max(0.25, 0.03*max(diff(x_limits), 1));
xlim(ax, x_limits + [-x_margin, x_margin])
xlabel(ax, "Feedback delay [ms]")
ylabel(ax, y_label)
title(ax, strrep(test.id, "_", " "))

if show_legend
    plot(ax, NaN, NaN, "kx", "LineWidth", 1.3, ...
        "DisplayName", "Saturation > 5%");
    plot(ax, NaN, NaN, "r*", "LineWidth", 1.1, ...
        "DisplayName", "Lost control");
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

