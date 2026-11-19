%% CALIBRATE_COMET_CYCLIC_EFFECTIVENESS
% One-parameter calibration of the global cyclic-effectiveness scale K_cyc.
%
% This script sweeps
%       P.aerodynamics.cyclic_effectiveness_scale
% while keeping the vehicle, controller, actuator limits and trajectory
% identical for the single- and dual-swashplate cases. Every candidate is
% repeated with one and two RK4 plant substeps while the controller remains
% fixed at 200 Hz, so numerical convergence is tested independently.
%
% The circle_4 experiment is used because it places the single configuration
% close to its cyclic-authority limit while the dual configuration remains
% controllable. This makes it the most informative case for identifying the
% static command-to-wrench effectiveness. Servo dynamics, delays, drag and
% disturbances must be calibrated later and are not hidden inside K_cyc.

clearvars
clc
close all

%% User settings

% The wide sweep identified the authority transition near 0.76--0.82. This
% local sweep resolves that interval while testing every candidate with one
% and two RK4 substeps at the same 200 Hz controller sample rate.
K_candidates = (0.760:0.005:0.820).';

calibration_case_id = "circle_4";
save_results = true;

% The score mainly identifies physical actuator demand. Tracking error has a
% deliberately small weight because the present nominal model does not yet
% contain the delays, uncertainty and disturbances seen in the experiment.
weights.command = 0.65;
weights.power = 0.25;
weights.tracking = 0.10;
convergence_tolerance_percent = 5.0;

%% Load the common model and experimental target

script_path = mfilename("fullpath");
script_folder = fileparts(script_path);
addpath(script_folder);

P_base = comet_parameters();
P_base.disturbance.enabled = false;
P_base.disturbance.type = "none";

assert(P_base.controller.same_gains_for_all_configurations, ...
    "The calibration requires identical controller gains.");
assert(all(isfinite(K_candidates)) && all(K_candidates > 0) && ...
       all(K_candidates <= 1), ...
    "Every K_cyc candidate must lie in the interval (0,1].");

test = select_case(P_base.experimental.tracking, calibration_case_id);

scenario = struct();
scenario.type = test.trajectory;
scenario.vmax = test.vmax;
scenario.amax = test.amax;
scenario.cycles = P_base.trajectory.default_cycles;
scenario.start_from_rest = true;

%% Preallocate scalar results

nK = numel(K_candidates);

single_rmse = nan(nK,1);
dual_rmse = nan(nK,1);
single_max_error = nan(nK,1);
dual_max_error = nan(nK,1);
single_velocity_rmse = nan(nK,1);
dual_velocity_rmse = nan(nK,1);

single_saturation = nan(nK,1);
dual_saturation = nan(nK,1);
single_position_saturation = nan(nK,1);
dual_position_saturation = nan(nK,1);
single_rate_saturation = nan(nK,1);
dual_rate_saturation = nan(nK,1);

single_peak_cyclic = nan(nK,1);
dual_peak_cyclic = nan(nK,1);
single_average_power = nan(nK,1);
dual_average_power = nan(nK,1);
single_maximum_power = nan(nK,1);
dual_maximum_power = nan(nK,1);

command_score = nan(nK,1);
power_score = nan(nK,1);
tracking_score = nan(nK,1);
total_score = nan(nK,1);

single_lost_control = false(nK,1);
dual_lost_control = false(nK,1);
single_lost_control_nominal = false(nK,1);
dual_lost_control_nominal = false(nK,1);

single_rmse_nominal = nan(nK,1);
dual_rmse_nominal = nan(nK,1);
single_max_error_nominal = nan(nK,1);
dual_max_error_nominal = nan(nK,1);
single_rmse_convergence_percent = nan(nK,1);
dual_rmse_convergence_percent = nan(nK,1);
single_max_error_convergence_percent = nan(nK,1);
dual_max_error_convergence_percent = nan(nK,1);
maximum_convergence_difference_percent = nan(nK,1);

fprintf("\nCOMET CYCLIC-EFFECTIVENESS CALIBRATION\n");
fprintf("Case: %s | %s | vmax = %.2f m/s\n", ...
    char(test.id), char(test.trajectory), test.vmax);
fprintf(['Only K_cyc and the internal RK4 substep count change; ' ...
         'the controller remains fixed at 200 Hz.\n\n']);

%% Sweep K_cyc

for i = 1:nK
    fprintf("[%2d/%2d] K_cyc = %.3f\n", i, nK, K_candidates(i));

    % Nominal integration: one RK4 step per 5 ms controller interval.
    P_nominal = P_base;
    P_nominal.aerodynamics.cyclic_effectiveness_scale = K_candidates(i);
    P_nominal.simulation.integration_substeps = 1;
    single_nominal = simulate_COMET(P_nominal, "single", scenario);
    dual_nominal = simulate_COMET(P_nominal, "dual", scenario);

    single_rmse_nominal(i) = single_nominal.metrics.rmse_position;
    dual_rmse_nominal(i) = dual_nominal.metrics.rmse_position;
    single_max_error_nominal(i) = ...
        single_nominal.metrics.max_position_error;
    dual_max_error_nominal(i) = ...
        dual_nominal.metrics.max_position_error;
    single_lost_control_nominal(i) = single_nominal.metrics.lost_control;
    dual_lost_control_nominal(i) = dual_nominal.metrics.lost_control;

    % Refined integration: two 2.5 ms RK4 steps under the same zero-order-
    % held command. Controller sampling, filters and integral updates remain
    % exactly at 200 Hz.
    P = P_base;
    P.aerodynamics.cyclic_effectiveness_scale = K_candidates(i);
    P.simulation.integration_substeps = 2;
    single = simulate_COMET(P, "single", scenario);
    dual = simulate_COMET(P, "dual", scenario);

    ms = single.metrics;
    md = dual.metrics;

    single_rmse(i) = ms.rmse_position;
    dual_rmse(i) = md.rmse_position;
    single_max_error(i) = ms.max_position_error;
    dual_max_error(i) = md.max_position_error;
    single_velocity_rmse(i) = ms.rmse_velocity;
    dual_velocity_rmse(i) = md.rmse_velocity;

    single_saturation(i) = 100*ms.any_saturation_fraction;
    dual_saturation(i) = 100*md.any_saturation_fraction;
    single_position_saturation(i) = ...
        100*ms.cyclic_position_saturation_fraction;
    dual_position_saturation(i) = ...
        100*md.cyclic_position_saturation_fraction;
    single_rate_saturation(i) = ...
        100*ms.actuator_rate_saturation_fraction;
    dual_rate_saturation(i) = ...
        100*md.actuator_rate_saturation_fraction;

    single_peak_cyclic(i) = peak_cyclic_use(single);
    dual_peak_cyclic(i) = peak_cyclic_use(dual);
    single_average_power(i) = ms.average_power_W;
    dual_average_power(i) = md.average_power_W;
    single_maximum_power(i) = ms.maximum_power_W;
    dual_maximum_power(i) = md.maximum_power_W;
    single_lost_control(i) = ms.lost_control;
    dual_lost_control(i) = md.lost_control;

    single_rmse_convergence_percent(i) = relative_difference_percent( ...
        single_rmse_nominal(i), single_rmse(i));
    dual_rmse_convergence_percent(i) = relative_difference_percent( ...
        dual_rmse_nominal(i), dual_rmse(i));
    single_max_error_convergence_percent(i) = ...
        relative_difference_percent( ...
        single_max_error_nominal(i), single_max_error(i));
    dual_max_error_convergence_percent(i) = ...
        relative_difference_percent( ...
        dual_max_error_nominal(i), dual_max_error(i));
    maximum_convergence_difference_percent(i) = max([ ...
        single_rmse_convergence_percent(i), ...
        dual_rmse_convergence_percent(i), ...
        single_max_error_convergence_percent(i), ...
        dual_max_error_convergence_percent(i)]);

    % Static cyclic demand is the primary evidence for K_cyc. Both ends of
    % the experimental roll/pitch command ranges are retained; using only an
    % absolute peak would discard directional asymmetry in the manoeuvre.
    simulated_commands = [ ...
        ms.roll_command_range, ms.pitch_command_range, ...
        md.roll_command_range, md.pitch_command_range];
    experimental_commands = [ ...
        test.single.roll_command_range, test.single.pitch_command_range, ...
        test.dual.roll_command_range, test.dual.pitch_command_range];
    command_score(i) = normalized_rms_mismatch( ...
        simulated_commands, experimental_commands, 0.10);

    simulated_power = [ ...
        ms.average_power_W, ms.maximum_power_W, ...
        md.average_power_W, md.maximum_power_W];
    experimental_power = [ ...
        test.single.average_power, test.single.maximum_power, ...
        test.dual.average_power, test.dual.maximum_power];
    power_score(i) = normalized_rms_mismatch( ...
        simulated_power, experimental_power, 50.0);

    simulated_tracking = [ ...
        ms.rmse_position, ms.max_position_error, ...
        md.rmse_position, md.max_position_error];
    experimental_tracking = [ ...
        test.single.rmse_position, test.single.max_position_error, ...
        test.dual.rmse_position, test.dual.max_position_error];
    tracking_score(i) = normalized_rms_mismatch( ...
        simulated_tracking, experimental_tracking, 0.05);

    total_score(i) = weights.command*command_score(i) + ...
                     weights.power*power_score(i) + ...
                     weights.tracking*tracking_score(i);

    % The full time histories are intentionally not retained: a sweep should
    % remain memory efficient even when the candidate grid is refined.
    clear single dual single_nominal dual_nominal
end

%% Build and print the calibration table

% A candidate must complete both integrations, remain controlled in both
% configurations and change by no more than the selected tolerance when the
% RK4 step is halved without changing the controller sample time.
valid_candidate = isfinite(total_score) & ...
    ~single_lost_control_nominal & ~dual_lost_control_nominal & ...
    ~single_lost_control & ~dual_lost_control & ...
    maximum_convergence_difference_percent <= ...
    convergence_tolerance_percent;

valid_indices = find(valid_candidate);
if isempty(valid_indices)
    error('calibrate_COMET:NoValidCandidate', ...
        ['No K_cyc candidate is simultaneously controlled and numerically ' ...
         'converged. Extend the local interval or use additional RK4 ' ...
         'substeps before calibrating the model.']);
end

calibration_table = table( ...
    K_candidates, ...
    single_rmse_nominal, single_rmse, ...
    dual_rmse_nominal, dual_rmse, ...
    single_max_error_nominal, single_max_error, ...
    dual_max_error_nominal, dual_max_error, ...
    single_rmse_convergence_percent, ...
    dual_rmse_convergence_percent, ...
    single_max_error_convergence_percent, ...
    dual_max_error_convergence_percent, ...
    maximum_convergence_difference_percent, ...
    single_velocity_rmse, dual_velocity_rmse, ...
    single_saturation, dual_saturation, ...
    single_position_saturation, dual_position_saturation, ...
    single_rate_saturation, dual_rate_saturation, ...
    single_peak_cyclic, dual_peak_cyclic, ...
    single_average_power, dual_average_power, ...
    single_maximum_power, dual_maximum_power, ...
    command_score, power_score, tracking_score, total_score, ...
    single_lost_control_nominal, dual_lost_control_nominal, ...
    single_lost_control, dual_lost_control, valid_candidate, ...
    'VariableNames', { ...
    'K_cyc', ...
    'Single_RMSE_1step_m', 'Single_RMSE_2step_m', ...
    'Dual_RMSE_1step_m', 'Dual_RMSE_2step_m', ...
    'Single_max_error_1step_m', 'Single_max_error_2step_m', ...
    'Dual_max_error_1step_m', 'Dual_max_error_2step_m', ...
    'Single_RMSE_difference_percent', ...
    'Dual_RMSE_difference_percent', ...
    'Single_max_error_difference_percent', ...
    'Dual_max_error_difference_percent', ...
    'Maximum_convergence_difference_percent', ...
    'Single_velocity_RMSE_mps', 'Dual_velocity_RMSE_mps', ...
    'Single_any_saturation_percent', 'Dual_any_saturation_percent', ...
    'Single_position_saturation_percent', ...
    'Dual_position_saturation_percent', ...
    'Single_rate_saturation_percent', 'Dual_rate_saturation_percent', ...
    'Single_peak_cyclic_percent', 'Dual_peak_cyclic_percent', ...
    'Single_average_power_W', 'Dual_average_power_W', ...
    'Single_maximum_power_W', 'Dual_maximum_power_W', ...
    'Command_score', 'Power_score', 'Tracking_score', 'Total_score', ...
    'Single_lost_control_1step', 'Dual_lost_control_1step', ...
    'Single_lost_control_2step', 'Dual_lost_control_2step', ...
    'Valid_candidate'});

% For an authority-bound calibration, select the lowest K_cyc that is both
% controlled and numerically converged. This preserves the experimentally
% relevant single-actuator saturation without selecting an unstable point.
recommended_index = valid_indices(1);

[~, local_index] = min(total_score(valid_candidate));
best_total_index = valid_indices(local_index);
[~, local_index] = min(command_score(valid_candidate));
best_command_index = valid_indices(local_index);
[~, local_index] = min(tracking_score(valid_candidate));
best_tracking_index = valid_indices(local_index);

fprintf("\nCALIBRATION RESULTS\n");
disp(calibration_table)

fprintf("Recommended K_cyc from the authority boundary: %.3f\n", ...
    K_candidates(recommended_index));
fprintf("Best K_cyc from the weighted physical score: %.3f\n", ...
    K_candidates(best_total_index));
fprintf("Best K_cyc from cyclic-command agreement only: %.3f\n", ...
    K_candidates(best_command_index));
fprintf("Best K_cyc from tracking agreement only: %.3f\n", ...
    K_candidates(best_tracking_index));

if recommended_index == 1 || recommended_index == nK
    warning('calibrate_COMET:BoundaryOptimum', ...
        ['The recommended value lies at the local sweep boundary. ' ...
         'Extend K_candidates before accepting the calibration.']);
end

fprintf(['\nDo not change comet_parameters.m automatically. ' ...
         'First inspect convergence and the plots; the provisional value is\n' ...
         'P.aerodynamics.cyclic_effectiveness_scale = %.3f;\n'], ...
    K_candidates(recommended_index));

%% Plot the parameter sweep

plot_score_sweep(K_candidates, command_score, power_score, ...
    tracking_score, total_score, K_candidates(recommended_index));

plot_convergence_sweep(K_candidates, ...
    single_rmse_convergence_percent, dual_rmse_convergence_percent, ...
    single_max_error_convergence_percent, ...
    dual_max_error_convergence_percent, ...
    convergence_tolerance_percent, K_candidates(recommended_index));

plot_tracking_sweep(K_candidates, single_rmse, dual_rmse, ...
    single_max_error, dual_max_error, test);

plot_actuation_sweep(K_candidates, single_peak_cyclic, dual_peak_cyclic, ...
    single_saturation, dual_saturation);

plot_power_sweep(K_candidates, single_average_power, dual_average_power, ...
    single_maximum_power, dual_maximum_power, test);

%% Save compact results

calibration = struct();
calibration.description = ...
    "Static cyclic-effectiveness sweep using the circle_4 experiment.";
calibration.case = test;
calibration.scenario = scenario;
calibration.candidates = K_candidates;
calibration.weights = weights;
calibration.table = calibration_table;
calibration.best_weighted_K = K_candidates(best_total_index);
calibration.best_command_K = K_candidates(best_command_index);
calibration.best_tracking_K = K_candidates(best_tracking_index);
calibration.recommended_K = K_candidates(recommended_index);
calibration.valid_candidates = valid_candidate;
calibration.convergence_tolerance_percent = ...
    convergence_tolerance_percent;
calibration.integration_substeps_compared = [1, 2];
calibration.base_parameters = P_base;

if save_results
    results_path = fullfile(script_folder, ...
        "COMET_cyclic_effectiveness_sweep.mat");
    save(results_path, "calibration");
    fprintf("\nCompact sweep results saved in:\n%s\n", results_path);
end

fprintf("\nCyclic-effectiveness sweep completed.\n");


%% Local functions

function test = select_case(cases, requested_id)

ids = string({cases.id});
index = find(ids == string(requested_id), 1);
if isempty(index)
    error('calibrate_COMET:UnknownCase', ...
        "Unknown case '%s'. Available cases: %s.", ...
        requested_id, strjoin(ids, ", "));
end
test = cases(index);

end


function score = normalized_rms_mismatch(simulated, experimental, floor_value)
% Dimensionless RMS mismatch. Missing paper entries are excluded.

simulated = simulated(:);
experimental = experimental(:);
valid = isfinite(simulated) & isfinite(experimental);
if ~any(valid)
    score = NaN;
    return
end

scale = max(abs(experimental(valid)), floor_value);
relative_error = (simulated(valid) - experimental(valid))./scale;
score = sqrt(mean(relative_error.^2));

end


function value = relative_difference_percent(nominal, refined)
% Relative change caused only by halving the RK4 plant-integration step.

value = 100*abs(nominal - refined)/max(abs(refined), 1e-9);

end


function value = peak_cyclic_use(result)
% Maximum actual active cyclic displacement as a percentage of its limit.

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


function plot_score_sweep(K, command_score, power_score, tracking_score, ...
    total_score, best_K)

figure('Color', 'w', 'Name', 'COMET Kcyc calibration score');
plot(K, command_score, '-o', 'LineWidth', 1.4, ...
    'DisplayName', 'Cyclic-command mismatch');
hold on
plot(K, power_score, '-s', 'LineWidth', 1.4, ...
    'DisplayName', 'Power mismatch');
plot(K, tracking_score, '-^', 'LineWidth', 1.4, ...
    'DisplayName', 'Tracking mismatch');
plot(K, total_score, '-d', 'LineWidth', 2.0, ...
    'DisplayName', 'Weighted score');
xline(best_K, '--k', sprintf('Selected K = %.3f', best_K), ...
    'LabelVerticalAlignment', 'bottom');
grid on
xlabel('Cyclic-effectiveness scale K_{cyc}')
ylabel('Normalized mismatch [-]')
title('Cyclic-effectiveness calibration objective')
legend('Location', 'best')

end


function plot_convergence_sweep(K, single_rmse, dual_rmse, ...
    single_maximum, dual_maximum, tolerance, selected_K)

figure('Color', 'w', 'Name', 'COMET RK4 convergence sweep');
plot(K, single_rmse, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Single RMSE');
hold on
plot(K, single_maximum, '-s', 'LineWidth', 1.5, ...
    'DisplayName', 'Single maximum error');
plot(K, dual_rmse, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Dual RMSE');
plot(K, dual_maximum, '-s', 'LineWidth', 1.5, ...
    'DisplayName', 'Dual maximum error');
yline(tolerance, '--r', sprintf('Tolerance = %.1f%%', tolerance), ...
    'LineWidth', 1.2);
xline(selected_K, '--k', sprintf('Selected K = %.3f', selected_K), ...
    'LabelVerticalAlignment', 'bottom');
grid on
xlabel('Cyclic-effectiveness scale K_{cyc}')
ylabel('One-step versus two-step difference [%]')
title('RK4 convergence at a fixed 200 Hz controller rate')
legend('Location', 'best')

end


function plot_tracking_sweep(K, single_rmse, dual_rmse, ...
    single_max_error, dual_max_error, test)

figure('Color', 'w', 'Name', 'COMET Kcyc tracking sweep');
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(K, single_rmse, '-o', 'LineWidth', 1.5, 'DisplayName', 'Simulation single');
hold on
plot(K, dual_rmse, '-o', 'LineWidth', 1.5, 'DisplayName', 'Simulation dual');
yline(test.single.rmse_position, '--', 'Experiment single', 'LineWidth', 1.2);
yline(test.dual.rmse_position, '--', 'Experiment dual', 'LineWidth', 1.2);
grid on
xlabel('K_{cyc}')
ylabel('Position RMSE [m]')
title('Tracking RMSE')
legend('Location', 'best')

nexttile
plot(K, single_max_error, '-s', 'LineWidth', 1.5, ...
    'DisplayName', 'Simulation single');
hold on
plot(K, dual_max_error, '-s', 'LineWidth', 1.5, ...
    'DisplayName', 'Simulation dual');
yline(test.single.max_position_error, '--', 'Experiment single', ...
    'LineWidth', 1.2);
yline(test.dual.max_position_error, '--', 'Experiment dual', ...
    'LineWidth', 1.2);
grid on
xlabel('K_{cyc}')
ylabel('Maximum position error [m]')
title('Maximum tracking error')
legend('Location', 'best')

sgtitle('circle 4.0 m/s: simulation versus experiment')

end


function plot_actuation_sweep(K, single_peak, dual_peak, ...
    single_saturation, dual_saturation)

figure('Color', 'w', 'Name', 'COMET Kcyc actuator sweep');
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(K, single_peak, '-o', 'LineWidth', 1.5, 'DisplayName', 'Single');
hold on
plot(K, dual_peak, '-o', 'LineWidth', 1.5, 'DisplayName', 'Dual');
yline(100, '--k', 'Position limit', 'LineWidth', 1.2);
grid on
xlabel('K_{cyc}')
ylabel('Peak actual cyclic use [% of limit]')
title('Required cyclic authority')
legend('Location', 'best')

nexttile
plot(K, single_saturation, '-s', 'LineWidth', 1.5, ...
    'DisplayName', 'Single');
hold on
plot(K, dual_saturation, '-s', 'LineWidth', 1.5, ...
    'DisplayName', 'Dual');
grid on
xlabel('K_{cyc}')
ylabel('Time in position or rate saturation [%]')
title('Saturation exposure')
legend('Location', 'best')

sgtitle('Actuator response while only K_{cyc} changes')

end


function plot_power_sweep(K, single_average, dual_average, ...
    single_maximum, dual_maximum, test)

figure('Color', 'w', 'Name', 'COMET Kcyc power sweep');
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(K, single_average, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Simulation single');
hold on
plot(K, dual_average, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Simulation dual');
yline(test.single.average_power, '--', 'Experiment single', 'LineWidth', 1.2);
yline(test.dual.average_power, '--', 'Experiment dual', 'LineWidth', 1.2);
grid on
xlabel('K_{cyc}')
ylabel('Average power [W]')
title('Average power')
legend('Location', 'best')

nexttile
plot(K, single_maximum, '-s', 'LineWidth', 1.5, ...
    'DisplayName', 'Simulation single');
hold on
plot(K, dual_maximum, '-s', 'LineWidth', 1.5, ...
    'DisplayName', 'Simulation dual');
yline(test.single.maximum_power, '--', 'Experiment single', 'LineWidth', 1.2);
yline(test.dual.maximum_power, '--', 'Experiment dual', 'LineWidth', 1.2);
grid on
xlabel('K_{cyc}')
ylabel('Maximum power [W]')
title('Maximum power')
legend('Location', 'best')

sgtitle('Power is diagnostic, not the primary K_{cyc} target')

end
