%% RUN_COMET_COMPARISON
% Baseline comparison between the single- and dual-swashplate COMET models.
%
% Required files in the same folder:
%   1. comet_parameters.m
%   2. simulate_COMET.m
%
% This script deliberately uses:
%   - the same physical parameters;
%   - the same trajectory references;
%   - the same controller architecture and gains;
%   - the same actuator limits;
% for both configurations. The only architectural change is the cyclic mask
% selected inside simulate_COMET.m. Therefore, any difference in tracking or
% saturation is caused by the available actuation authority, not by a more
% favourable controller tuning for one configuration.
%
% The script performs the five trajectory cases reported in the paper,
% builds numerical comparison tables and plots:
%   - all ten simulated trajectories in one 2-by-5 window;
%   - grouped bar charts for tracking performance;
%   - grouped bar charts for actuator use and saturation;
%   - grouped bar charts comparing simulation and paper experiments.
%
% Position RMSE and maximum error are reported both physically and with the
% selected 60 ms measurement-timestamp lag. The lag is a post-processing
% model of acquisition/log synchronisation and is never fed to the controller.

clearvars
clc
close all

%% User settings

% Select the cases to simulate. Available identifiers are stored in
% P.experimental.tracking and are currently:
%   circle_3, circle_4, eight_3p5, eight_4, eight_5.
case_ids = ["circle_3", "circle_4", "eight_3p5", "eight_4", "eight_5"];

% Set true to save the complete MATLAB results structure beside this file.
save_results = true;

%% Load the common model

script_path = mfilename("fullpath");
script_folder = fileparts(script_path);
addpath(script_folder);

P = comet_parameters();

assert(abs(P.aerodynamics.cyclic_effectiveness_scale - 0.780) < 1e-12, ...
    "run_COMET_comparison expects K_cyc = 0.780.");
assert(abs(P.actuators.servo_time_constant - 0.020) < 1e-12, ...
    "run_COMET_comparison expects servo tau = 20 ms.");
assert(abs(P.feedback.state_delay - 0.0025) < 1e-12, ...
    "run_COMET_comparison expects a 2.5 ms physical feedback delay.");
assert(abs(P.trajectory.entry_time_override - 1.30) < 1e-12, ...
    "run_COMET_comparison expects a 1.30 s entry ramp.");
assert(isfield(P, "measurement") && ...
       isfield(P.measurement, "position_timestamp_lag") && ...
       abs(P.measurement.position_timestamp_lag - 0.060) < 1e-12, ...
    "run_COMET_comparison expects a 60 ms measurement timestamp lag.");

% This file is the undisturbed baseline. Disturbance tests will be performed
% in a separate file, after the nominal single-dual comparison is verified.
P.disturbance.enabled = false;
P.disturbance.type = "none";

assert(P.controller.same_gains_for_all_configurations, ...
    "Single and dual simulations must use the same controller gains.");
assert(P.allocation.configuration_dependent_cyclic_map, ...
    "The comparison requires configuration-dependent cyclic allocation.");

all_cases = P.experimental.tracking;
selected_cases = select_cases(all_cases, case_ids);
n_cases = numel(selected_cases);

%% Run the single and dual configurations

empty_run = struct("case_id", "", "scenario", struct(), ...
    "single", struct(), "dual", struct());
runs = repmat(empty_run, n_cases, 1);

fprintf("\nCOMET SINGLE-DUAL TRACKING COMPARISON\n");
fprintf("Same model, controller gains and limits for both configurations.\n\n");
fprintf(['K_cyc = %.3f | servo tau = %.1f ms | physical delay = %.1f ms | ' ...
    'measurement lag = %.0f ms | entry ramp = %.2f s\n\n'], ...
    P.aerodynamics.cyclic_effectiveness_scale, ...
    1000*P.actuators.servo_time_constant, ...
    1000*P.feedback.state_delay, ...
    1000*P.measurement.position_timestamp_lag, ...
    P.trajectory.entry_time_override);

for k = 1:n_cases
    test = selected_cases(k);

    scenario = struct();
    scenario.type = test.trajectory;
    scenario.vmax = test.vmax;
    scenario.amax = test.amax;
    scenario.cycles = P.trajectory.default_cycles;
    scenario.start_from_rest = P.trajectory.start_from_rest;
    scenario.initial_velocity_I = zeros(3,1);
    if isfinite(P.trajectory.entry_time_override)
        scenario.entry_time = P.trajectory.entry_time_override;
    end

    fprintf("[%d/%d] %-11s | %s | vmax = %.2f m/s\n", ...
        k, n_cases, char(test.id), char(test.trajectory), test.vmax);

    % Exactly the same P and scenario are passed to the two simulations.
    single_result = simulate_COMET(P, "single", scenario);
    dual_result = simulate_COMET(P, "dual", scenario);

    % Preserve the physical metrics and add the apparent metrics produced by
    % the selected measurement timestamp lag. Dynamics and saturation are
    % already complete and are not modified here.
    single_result = attach_reported_position_metrics( ...
        single_result, P.measurement.position_timestamp_lag);
    dual_result = attach_reported_position_metrics( ...
        dual_result, P.measurement.position_timestamp_lag);

    assert(abs(single_result.parameters.feedback.state_delay_applied - ...
               P.feedback.state_delay) < 1e-12 && ...
           abs(dual_result.parameters.feedback.state_delay_applied - ...
               P.feedback.state_delay) < 1e-12, ...
        "simulate_COMET did not preserve the selected physical feedback delay.");

    runs(k).case_id = test.id;
    runs(k).scenario = scenario;
    runs(k).single = single_result;
    runs(k).dual = dual_result;
end

delay_check = runs(1).single.parameters.feedback;
fprintf(['\nPhysical feedback-delay check: requested %.1f ms | ' ...
    'applied %.1f ms | floor samples %.0f | fraction %.3f\n'], ...
    1000*delay_check.state_delay_requested, ...
    1000*delay_check.state_delay_applied, ...
    delay_check.state_delay_floor_samples, ...
    delay_check.state_delay_fraction);

%% Build the numerical comparison tables

summary_table = build_summary_table(selected_cases, runs);
advantage_table = build_advantage_table(selected_cases, runs);
paper_table = build_paper_table(selected_cases, runs);

fprintf("\nSIMULATION SUMMARY\n");
disp(summary_table)

fprintf("\nDUAL CONFIGURATION ADVANTAGE\n");
disp(advantage_table)

single_saturation_cases = ...
    advantage_table.Case(advantage_table.Single_saturation_detected);
if isempty(single_saturation_cases)
    fprintf("No single-configuration actuator saturation detected.\n");
else
    fprintf("Single-configuration saturation detected in: %s\n", ...
        strjoin(single_saturation_cases, ", "));
end

fprintf("\nSIMULATION VERSUS REPORTED EXPERIMENTAL RESULTS\n");
disp(paper_table)

%% Plot all trajectories and categorical comparisons

plot_all_trajectories(selected_cases, runs);
plot_tracking_bars(selected_cases, runs);
plot_actuation_bars(selected_cases, runs);
plot_paper_bars(selected_cases, runs);

%% Collect and optionally save all results

study = struct();
study.description = strcat( ...
    "Nominal COMET single-dual comparison with identical control settings ", ...
    "and a post-processed position measurement timestamp lag.");
study.parameters = P;
study.measurement_timestamp_lag_s = ...
    P.measurement.position_timestamp_lag;
study.cases = selected_cases;
study.runs = runs;
study.summary = summary_table;
study.dual_advantage = advantage_table;
study.paper_comparison = paper_table;

if save_results
    results_path = fullfile(script_folder, "COMET_comparison_results.mat");
    save(results_path, "study", "-v7.3");
    fprintf("\nComplete results saved in:\n%s\n", results_path);
end

fprintf("\nComparison completed.\n");


%% Local functions

function selected = select_cases(all_cases, requested_ids)
% Return the requested experimental cases in the requested order.

requested_ids = string(requested_ids(:));
available_ids = string([all_cases.id]).';

selected = repmat(all_cases(1), numel(requested_ids), 1);
for i = 1:numel(requested_ids)
    index = find(available_ids == requested_ids(i), 1);
    if isempty(index)
        error("run_COMET_comparison:UnknownCase", ...
            "Unknown case '%s'. Available cases: %s.", ...
            requested_ids(i), strjoin(available_ids, ", "));
    end
    selected(i) = all_cases(index);
end

end


function T = build_summary_table(cases, runs)
% One row per case and configuration.

n = 2*numel(cases);
Case = strings(n,1);
Trajectory = strings(n,1);
Configuration = strings(n,1);
Vmax_mps = nan(n,1);
Amax_mps2 = nan(n,1);
Measurement_lag_ms = nan(n,1);
RMSE_position_reported_m = nan(n,1);
RMSE_position_physical_m = nan(n,1);
Max_position_error_reported_m = nan(n,1);
Max_position_error_physical_m = nan(n,1);
RMSE_velocity_mps = nan(n,1);
Any_saturation_percent = nan(n,1);
Saturation_detected = false(n,1);
Cyclic_position_sat_percent = nan(n,1);
Rate_saturation_percent = nan(n,1);
Peak_cyclic_use_percent = nan(n,1);
Minimum_allocation_rcond = nan(n,1);
Allocation_pseudoinverse_percent = nan(n,1);
Average_power_W = nan(n,1);
Maximum_power_W = nan(n,1);
Lost_control = false(n,1);

row = 0;
for k = 1:numel(cases)
    names = ["single", "dual"];
    for j = 1:2
        row = row + 1;
        result = runs(k).(names(j));
        m = result.metrics;

        Case(row) = cases(k).id;
        Trajectory(row) = cases(k).trajectory;
        Configuration(row) = names(j);
        Vmax_mps(row) = cases(k).vmax;
        Amax_mps2(row) = cases(k).amax;
        Measurement_lag_ms(row) = 1000*m.position_timestamp_lag_s;
        RMSE_position_reported_m(row) = m.rmse_position_reported;
        RMSE_position_physical_m(row) = m.rmse_position_physical;
        Max_position_error_reported_m(row) = ...
            m.max_position_error_reported;
        Max_position_error_physical_m(row) = ...
            m.max_position_error_physical;
        RMSE_velocity_mps(row) = m.rmse_velocity;
        Any_saturation_percent(row) = 100*m.any_saturation_fraction;
        Saturation_detected(row) = m.any_saturation_fraction > 0;
        Cyclic_position_sat_percent(row) = ...
            100*m.cyclic_position_saturation_fraction;
        Rate_saturation_percent(row) = ...
            100*m.actuator_rate_saturation_fraction;
        Peak_cyclic_use_percent(row) = peak_cyclic_use(result);
        Minimum_allocation_rcond(row) = ...
            m.minimum_allocation_reciprocal_condition;
        Allocation_pseudoinverse_percent(row) = ...
            100*m.allocation_pseudoinverse_fraction;
        Average_power_W(row) = m.average_power_W;
        Maximum_power_W(row) = m.maximum_power_W;
        Lost_control(row) = m.lost_control;
    end
end

T = table(Case, Trajectory, Configuration, Vmax_mps, Amax_mps2, ...
    Measurement_lag_ms, ...
    RMSE_position_reported_m, RMSE_position_physical_m, ...
    Max_position_error_reported_m, Max_position_error_physical_m, ...
    RMSE_velocity_mps, ...
    Any_saturation_percent, Saturation_detected, ...
    Cyclic_position_sat_percent, ...
    Rate_saturation_percent, Peak_cyclic_use_percent, ...
    Minimum_allocation_rcond, Allocation_pseudoinverse_percent, ...
    Average_power_W, Maximum_power_W, Lost_control);

end


function T = build_advantage_table(cases, runs)
% Positive reductions mean that the dual configuration performs better.

n = numel(cases);
Case = strings(n,1);
Trajectory = strings(n,1);
Vmax_mps = nan(n,1);
Reported_RMSE_reduction_percent = nan(n,1);
Physical_RMSE_reduction_percent = nan(n,1);
Reported_max_error_reduction_percent = nan(n,1);
Physical_max_error_reduction_percent = nan(n,1);
Velocity_RMSE_reduction_percent = nan(n,1);
Single_saturation_percent = nan(n,1);
Dual_saturation_percent = nan(n,1);
Single_saturation_detected = false(n,1);
Dual_saturation_detected = false(n,1);
Single_peak_cyclic_percent = nan(n,1);
Dual_peak_cyclic_percent = nan(n,1);
Single_lost_control = false(n,1);
Dual_lost_control = false(n,1);

for k = 1:n
    s = runs(k).single.metrics;
    d = runs(k).dual.metrics;

    Case(k) = cases(k).id;
    Trajectory(k) = cases(k).trajectory;
    Vmax_mps(k) = cases(k).vmax;
    Reported_RMSE_reduction_percent(k) = percent_reduction( ...
        s.rmse_position_reported, d.rmse_position_reported);
    Physical_RMSE_reduction_percent(k) = percent_reduction( ...
        s.rmse_position_physical, d.rmse_position_physical);
    Reported_max_error_reduction_percent(k) = percent_reduction( ...
        s.max_position_error_reported, d.max_position_error_reported);
    Physical_max_error_reduction_percent(k) = percent_reduction( ...
        s.max_position_error_physical, d.max_position_error_physical);
    Velocity_RMSE_reduction_percent(k) = percent_reduction( ...
        s.rmse_velocity, d.rmse_velocity);
    Single_saturation_percent(k) = 100*s.any_saturation_fraction;
    Dual_saturation_percent(k) = 100*d.any_saturation_fraction;
    Single_saturation_detected(k) = s.any_saturation_fraction > 0;
    Dual_saturation_detected(k) = d.any_saturation_fraction > 0;
    Single_peak_cyclic_percent(k) = peak_cyclic_use(runs(k).single);
    Dual_peak_cyclic_percent(k) = peak_cyclic_use(runs(k).dual);
    Single_lost_control(k) = s.lost_control;
    Dual_lost_control(k) = d.lost_control;
end

T = table(Case, Trajectory, Vmax_mps, ...
    Reported_RMSE_reduction_percent, Physical_RMSE_reduction_percent, ...
    Reported_max_error_reduction_percent, ...
    Physical_max_error_reduction_percent, ...
    Velocity_RMSE_reduction_percent, ...
    Single_saturation_percent, Dual_saturation_percent, ...
    Single_saturation_detected, Dual_saturation_detected, ...
    Single_peak_cyclic_percent, Dual_peak_cyclic_percent, ...
    Single_lost_control, Dual_lost_control);

end


function T = build_paper_table(cases, runs)
% Compare the simulation with values explicitly reported in the paper.

n = 2*numel(cases);
Case = strings(n,1);
Configuration = strings(n,1);
Experimental_available = false(n,1);
Experimental_completed = false(n,1);
Measurement_lag_ms = nan(n,1);
Simulation_reported_RMSE_m = nan(n,1);
Simulation_physical_RMSE_m = nan(n,1);
Experimental_RMSE_m = nan(n,1);
RMSE_difference_m = nan(n,1);
Simulation_reported_max_error_m = nan(n,1);
Simulation_physical_max_error_m = nan(n,1);
Experimental_max_error_m = nan(n,1);
Max_error_difference_m = nan(n,1);

row = 0;
for k = 1:numel(cases)
    names = ["single", "dual"];
    for j = 1:2
        row = row + 1;
        name = names(j);
        target = cases(k).(name);
        sim = runs(k).(name).metrics;

        Case(row) = cases(k).id;
        Configuration(row) = name;
        Experimental_available(row) = target.available;
        Experimental_completed(row) = target.completed;
        Measurement_lag_ms(row) = 1000*sim.position_timestamp_lag_s;
        Simulation_reported_RMSE_m(row) = sim.rmse_position_reported;
        Simulation_physical_RMSE_m(row) = sim.rmse_position_physical;
        Experimental_RMSE_m(row) = target.rmse_position;
        Simulation_reported_max_error_m(row) = ...
            sim.max_position_error_reported;
        Simulation_physical_max_error_m(row) = ...
            sim.max_position_error_physical;
        Experimental_max_error_m(row) = target.max_position_error;

        if target.available && target.completed
            RMSE_difference_m(row) = ...
                sim.rmse_position_reported - target.rmse_position;
            Max_error_difference_m(row) = ...
                sim.max_position_error_reported - ...
                target.max_position_error;
        end
    end
end

T = table(Case, Configuration, Experimental_available, ...
    Experimental_completed, Measurement_lag_ms, ...
    Simulation_reported_RMSE_m, Simulation_physical_RMSE_m, ...
    Experimental_RMSE_m, RMSE_difference_m, ...
    Simulation_reported_max_error_m, Simulation_physical_max_error_m, ...
    Experimental_max_error_m, Max_error_difference_m);

end


function result = attach_reported_position_metrics(result, timestamp_lag_s)
% Add post-processed tracking metrics without changing physical dynamics.

[reported_rmse, reported_max_error] = ...
    apparent_position_metrics(result, timestamp_lag_s);

result.metrics.position_timestamp_lag_s = timestamp_lag_s;
result.metrics.rmse_position_physical = result.metrics.rmse_position;
result.metrics.max_position_error_physical = ...
    result.metrics.max_position_error;
result.metrics.rmse_position_reported = reported_rmse;
result.metrics.max_position_error_reported = reported_max_error;

end


function [rmse_position, max_position_error] = ...
    apparent_position_metrics(result, timestamp_lag_s)
% Compare r(t) with the position measurement x(t-lag).

t = result.time(:);
query_time = max(t - timestamp_lag_s, t(1));
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


function value = percent_reduction(single_value, dual_value)

if ~isfinite(single_value) || abs(single_value) < eps
    value = NaN;
else
    value = 100*(single_value - dual_value)/single_value;
end

end


function value = peak_cyclic_use(result)
% Maximum active cyclic displacement as a percentage of its hard limit.

evaluation = result.time >= result.scenario.metric_start_time;
active_rows = 2 + find(logical(result.cyclic_mask(:)));
cyclic = result.actuator.actual(active_rows, evaluation);
limit = max(abs([result.parameters.actuators.cyclic_min, ...
                 result.parameters.actuators.cyclic_max]));

if isempty(cyclic) || limit <= 0
    value = 0;
else
    value = 100*max(abs(cyclic), [], "all")/limit;
end

end


function labels = display_case_labels(cases)
% Readable tick labels without TeX subscripts.

labels = string([cases.id]);
labels = strrep(labels, "_3p5", " 3.5");
labels = strrep(labels, "_", " ");

end


function plot_all_trajectories(cases, runs)
% Show every simulated single and dual trajectory in one figure.

n = numel(cases);
single_color = [0.0000, 0.4470, 0.7410];
dual_color = [0.8500, 0.3250, 0.0980];

figure("Name", "COMET all simulated trajectories", "Color", "w", ...
    "Position", [60, 80, 1750, 720]);
tiledlayout(2,n, "TileSpacing", "compact", "Padding", "compact");

for row = 1:2
    for k = 1:n
        if row == 1
            result = runs(k).single;
            configuration = "Single";
            actual_color = single_color;
        else
            result = runs(k).dual;
            configuration = "Dual";
            actual_color = dual_color;
        end

        nexttile
        plot(result.reference.position_I(1,:), ...
            result.reference.position_I(2,:), "k--", "LineWidth", 1.3)
        hold on
        plot(result.state.position_I(1,:), result.state.position_I(2,:), ...
            "Color", actual_color, "LineWidth", 1.5)
        axis equal
        grid on
        xlabel("x_I [m]")
        if k == 1
            ylabel(configuration + "   y_I [m]")
        end
        title(sprintf("%s | %.1f m/s", ...
            char(cases(k).trajectory), cases(k).vmax))

        if k == 1
            legend("Reference", configuration, "Location", "best")
        end
    end
end

sgtitle("All nominal COMET trajectory-tracking simulations")

end


function plot_tracking_bars(cases, runs)
% Categorical tracking metrics: grouped bars are more appropriate than lines.

n = numel(cases);
x = 1:n;
labels = display_case_labels(cases);

single_rmse = arrayfun(@(r) ...
    r.single.metrics.rmse_position_reported, runs);
dual_rmse = arrayfun(@(r) ...
    r.dual.metrics.rmse_position_reported, runs);
single_max = arrayfun(@(r) ...
    r.single.metrics.max_position_error_reported, runs);
dual_max = arrayfun(@(r) ...
    r.dual.metrics.max_position_error_reported, runs);
single_velocity = arrayfun(@(r) r.single.metrics.rmse_velocity, runs);
dual_velocity = arrayfun(@(r) r.dual.metrics.rmse_velocity, runs);

figure("Name", "COMET tracking performance bars", "Color", "w", ...
    "Position", [90, 120, 1550, 470]);
tiledlayout(1,3, "TileSpacing", "compact", "Padding", "compact");

nexttile
b = bar(x, [single_rmse(:), dual_rmse(:)], "grouped");
set_pair_colors(b)
grid on
set(gca, "XTick", x, "XTickLabel", labels)
ylabel("Position RMSE [m]")
title("Reported position RMSE (60 ms lag)")
legend("Single", "Dual", "Location", "northwest")

nexttile
b = bar(x, [single_max(:), dual_max(:)], "grouped");
set_pair_colors(b)
grid on
set(gca, "XTick", x, "XTickLabel", labels)
ylabel("Maximum position error [m]")
title("Reported maximum error (60 ms lag)")

nexttile
b = bar(x, [single_velocity(:), dual_velocity(:)], "grouped");
set_pair_colors(b)
grid on
set(gca, "XTick", x, "XTickLabel", labels)
ylabel("Velocity RMSE [m/s]")
title("Velocity RMSE")

sgtitle(['Tracking performance with identical controller gains; ' ...
    'position metrics include the selected timestamp lag'])

end


function plot_actuation_bars(cases, runs)
% Compare demanded cyclic authority and actual saturation occurrence.

n = numel(cases);
x = 1:n;
labels = display_case_labels(cases);

single_peak = arrayfun(@(r) peak_cyclic_use(r.single), runs);
dual_peak = arrayfun(@(r) peak_cyclic_use(r.dual), runs);
single_sat = 100*arrayfun(@(r) ...
    r.single.metrics.any_saturation_fraction, runs);
dual_sat = 100*arrayfun(@(r) ...
    r.dual.metrics.any_saturation_fraction, runs);

figure("Name", "COMET actuator use bars", "Color", "w", ...
    "Position", [140, 150, 1200, 470]);
tiledlayout(1,2, "TileSpacing", "compact", "Padding", "compact");

nexttile
b = bar(x, [single_peak(:), dual_peak(:)], "grouped");
set_pair_colors(b)
hold on
yline(100, "k--", "Actuator limit", "LineWidth", 1.2)
grid on
set(gca, "XTick", x, "XTickLabel", labels)
ylabel("Peak cyclic use [% of limit]")
title("Maximum cyclic demand")
legend("Single", "Dual", "Limit", "Location", "northwest")

nexttile
b = bar(x, [single_sat(:), dual_sat(:)], "grouped");
set_pair_colors(b)
grid on
set(gca, "XTick", x, "XTickLabel", labels)
ylabel("Time in saturation [%]")
title("Position or rate saturation")

if all([single_sat(:); dual_sat(:)] == 0)
    ylim([0, 1])
    text(mean(xlim), 0.55, ...
        "No actuator saturation in any nominal case", ...
        "HorizontalAlignment", "center", "FontWeight", "bold")
end

sgtitle("Actuator use with identical limits")

end


function plot_paper_bars(cases, runs)
% Compare every reported paper metric that has a simulated counterpart.

n = numel(cases);
x = 1:n;
labels = display_case_labels(cases);

sim_rmse_single = arrayfun(@(r) ...
    r.single.metrics.rmse_position_reported, runs);
sim_rmse_dual = arrayfun(@(r) ...
    r.dual.metrics.rmse_position_reported, runs);
sim_max_single = arrayfun(@(r) ...
    r.single.metrics.max_position_error_reported, runs);
sim_max_dual = arrayfun(@(r) ...
    r.dual.metrics.max_position_error_reported, runs);
sim_pavg_single = arrayfun(@(r) r.single.metrics.average_power_W, runs);
sim_pavg_dual = arrayfun(@(r) r.dual.metrics.average_power_W, runs);
sim_pmax_single = arrayfun(@(r) r.single.metrics.maximum_power_W, runs);
sim_pmax_dual = arrayfun(@(r) r.dual.metrics.maximum_power_W, runs);
sim_roll_single = arrayfun(@(r) ...
    max(abs(r.single.metrics.roll_command_range)), runs);
sim_roll_dual = arrayfun(@(r) ...
    max(abs(r.dual.metrics.roll_command_range)), runs);
sim_pitch_single = arrayfun(@(r) ...
    max(abs(r.single.metrics.pitch_command_range)), runs);
sim_pitch_dual = arrayfun(@(r) ...
    max(abs(r.dual.metrics.pitch_command_range)), runs);

exp_rmse_single = nan(1,n);
exp_rmse_dual = nan(1,n);
exp_max_single = nan(1,n);
exp_max_dual = nan(1,n);
exp_pavg_single = nan(1,n);
exp_pavg_dual = nan(1,n);
exp_pmax_single = nan(1,n);
exp_pmax_dual = nan(1,n);
exp_roll_single = nan(1,n);
exp_roll_dual = nan(1,n);
exp_pitch_single = nan(1,n);
exp_pitch_dual = nan(1,n);

for k = 1:n
    if cases(k).single.available && cases(k).single.completed
        target = cases(k).single;
        exp_rmse_single(k) = target.rmse_position;
        exp_max_single(k) = target.max_position_error;
        exp_pavg_single(k) = target.average_power;
        exp_pmax_single(k) = target.maximum_power;
        exp_roll_single(k) = max(abs(target.roll_command_range));
        exp_pitch_single(k) = max(abs(target.pitch_command_range));
    end
    if cases(k).dual.available && cases(k).dual.completed
        target = cases(k).dual;
        exp_rmse_dual(k) = target.rmse_position;
        exp_max_dual(k) = target.max_position_error;
        exp_pavg_dual(k) = target.average_power;
        exp_pmax_dual(k) = target.maximum_power;
        exp_roll_dual(k) = max(abs(target.roll_command_range));
        exp_pitch_dual(k) = max(abs(target.pitch_command_range));
    end
end

figure("Name", "COMET simulation versus experiment bars", ...
    "Color", "w", "Position", [40, 50, 1750, 850]);
tiledlayout(2,3, "TileSpacing", "compact", "Padding", "compact");

nexttile
plot_experiment_bars(x, labels, ...
    [sim_rmse_single(:), exp_rmse_single(:), ...
     sim_rmse_dual(:), exp_rmse_dual(:)], ...
    "Position RMSE [m]", "Position RMSE", true)

nexttile
plot_experiment_bars(x, labels, ...
    [sim_max_single(:), exp_max_single(:), ...
     sim_max_dual(:), exp_max_dual(:)], ...
    "Maximum error [m]", "Maximum position error", false)

nexttile
plot_experiment_bars(x, labels, ...
    [sim_pavg_single(:), exp_pavg_single(:), ...
     sim_pavg_dual(:), exp_pavg_dual(:)], ...
    "Average power [W]", "Average power", false)

nexttile
plot_experiment_bars(x, labels, ...
    [sim_pmax_single(:), exp_pmax_single(:), ...
     sim_pmax_dual(:), exp_pmax_dual(:)], ...
    "Maximum power [W]", "Maximum power", false)

nexttile
plot_experiment_bars(x, labels, ...
    [sim_roll_single(:), exp_roll_single(:), ...
     sim_roll_dual(:), exp_roll_dual(:)], ...
    "Peak |roll command| [-]", "Roll-command demand", false)

nexttile
plot_experiment_bars(x, labels, ...
    [sim_pitch_single(:), exp_pitch_single(:), ...
     sim_pitch_dual(:), exp_pitch_dual(:)], ...
    "Peak |pitch command| [-]", "Pitch-command demand", false)

sgtitle(['Numerical reproduction versus reported experiments; ' ...
    'position metrics include the selected 60 ms timestamp lag'])

end


function plot_experiment_bars(x, labels, values, y_label, title_text, ...
    show_legend)
% Four bars per case: simulation/experiment for single/dual architectures.

b = bar(x, values, "grouped");
b(1).FaceColor = [0.0000, 0.4470, 0.7410];
b(2).FaceColor = [0.3010, 0.7450, 0.9330];
b(3).FaceColor = [0.8500, 0.3250, 0.0980];
b(4).FaceColor = [0.9290, 0.6940, 0.1250];
grid on
set(gca, "XTick", x, "XTickLabel", labels)
ylabel(y_label)
title(title_text)

if show_legend
    legend("Simulation single", "Experiment single", ...
        "Simulation dual", "Experiment dual", ...
        "Location", "northwest")
end

end


function set_pair_colors(bar_handles)
% Consistent colors in every single-dual grouped bar chart.

bar_handles(1).FaceColor = [0.0000, 0.4470, 0.7410];
bar_handles(2).FaceColor = [0.8500, 0.3250, 0.0980];

end
