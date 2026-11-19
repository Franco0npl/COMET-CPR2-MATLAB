%% CALIBRATE_COMET_MULTIPARAMETER
% Targeted two-parameter search for COMET tracking fidelity.
%
% Six physically plausible combinations of feedback delay and horizontal
% quadratic drag are tested. The cyclic-servo time constant is fixed at
% 30 ms because the previous one-parameter study showed little influence on
% stable cases and the least pathological circle-4-single response at 30 ms.
%
% Selection uses only:
%   1) position RMSE;
%   2) maximum position error.
%
% Only completed experiments with numerical values enter the score. Single
% eight 4 is checked separately as an expected loss-of-control case; single
% eight 5 is excluded because its experimental result was not reported.

clearvars
clc
close all

%% Targeted candidate definitions

fixed_servo_tau_s = 0.030;

candidate_labels = ["C1"; "C2"; "C3"; "C4"; "C5"; "C6"];
candidate_delay_s = [0.005; 0.005; 0.005; 0.010; 0.010; 0.010];
candidate_drag_xy = [0.0450; 0.0525; 0.0600; 0.0150; 0.0250; 0.0350];

acceptable_mismatch_percent = 10.0;
completed_case_saturation_guardrail_percent = 5.0;

save_results = true;
resume_from_checkpoint = true;

%% Load the fixed model

script_path = mfilename("fullpath");
script_folder = fileparts(script_path);
addpath(script_folder);

clear simulate_COMET comet_parameters
rehash

simulator_path = which("simulate_COMET");
parameter_path = which("comet_parameters");
fprintf("Using simulate_COMET: %s\n", simulator_path);
fprintf("Using comet_parameters: %s\n\n", parameter_path);

P_base = comet_parameters();
P_base.disturbance.enabled = false;
P_base.disturbance.type = "none";

assert(abs(P_base.aerodynamics.cyclic_effectiveness_scale - 0.78) < 1e-12, ...
    "This search expects K_cyc = 0.78 in comet_parameters.m.");
assert(P_base.controller.same_gains_for_all_configurations, ...
    "Single and dual must use identical controller gains.");

controller_dt = P_base.simulation.controller_sample_time;
candidate_delay_samples = round(candidate_delay_s/controller_dt);
quantized_delays_s = candidate_delay_samples*controller_dt;

assert(all(abs(quantized_delays_s - candidate_delay_s) < 1e-12), ...
    sprintf(["Every delay must be representable by the fixed controller " + ...
    "sample time (currently %.3f s)."], controller_dt));
candidate_delay_s = quantized_delays_s;

assert(fixed_servo_tau_s >= 0.030 && fixed_servo_tau_s <= 0.050, ...
    "The fixed servo time constant must remain inside 0.030--0.050 s.");
assert(all(candidate_delay_s >= 0.005 & candidate_delay_s <= 0.010), ...
    "Targeted delay values must remain inside 0.005--0.010 s.");
assert(all(candidate_drag_xy >= 0.015 & candidate_drag_xy <= 0.060), ...
    "Quadratic drag must remain inside 0.015--0.060 N/(m/s)^2.");
assert(numel(candidate_labels) == numel(candidate_delay_s) && ...
       numel(candidate_labels) == numel(candidate_drag_xy), ...
    "Candidate arrays must have identical lengths.");

n_candidates = numel(candidate_labels);
candidate_tau_s = fixed_servo_tau_s*ones(n_candidates,1);

candidates = table((1:n_candidates).', candidate_labels, ...
    1000*candidate_tau_s, 1000*candidate_delay_s, candidate_drag_xy, ...
    'VariableNames', {'Candidate', 'Label', 'Servo_tau_ms', ...
    'Feedback_delay_ms', 'Quadratic_drag_xy'});

tests = P_base.experimental.tracking;
runs = make_run_list(tests);
n_runs = numel(runs);
n_rows = n_candidates*n_runs;

fprintf("COMET TARGETED DELAY/DRAG COMBINATION SEARCH\n");
fprintf("K_cyc fixed at %.3f | servo tau fixed at %.0f ms\n", ...
    P_base.aerodynamics.cyclic_effectiveness_scale, ...
    1000*fixed_servo_tau_s);
disp(candidates)
fprintf("Candidates: %d | simulations per candidate: %d | total: %d\n\n", ...
    n_candidates, n_runs, n_rows);

%% Results folder and checkpoint

results_folder = fullfile(script_folder, "targeted_delay_drag_results");
if save_results && ~isfolder(results_folder)
    mkdir(results_folder)
end

checkpoint_path = fullfile(results_folder, ...
    "COMET_targeted_delay_drag_checkpoint.mat");
checkpoint_signature = make_checkpoint_signature(P_base, script_folder);

detailed_table = empty_detailed_table(n_rows);
candidate_complete = false(n_candidates,1);

if resume_from_checkpoint && isfile(checkpoint_path)
    saved = load(checkpoint_path);
    compatible = isfield(saved, "candidate_tau_s") && ...
        isfield(saved, "candidate_delay_s") && ...
        isfield(saved, "candidate_drag_xy") && ...
        isfield(saved, "candidate_labels") && ...
        isfield(saved, "checkpoint_signature") && ...
        isfield(saved, "detailed_table") && ...
        isfield(saved, "candidate_complete") && ...
        isequal(saved.candidate_tau_s, candidate_tau_s) && ...
        isequal(saved.candidate_delay_s, candidate_delay_s) && ...
        isequal(saved.candidate_drag_xy, candidate_drag_xy) && ...
        isequal(saved.candidate_labels, candidate_labels) && ...
        isequaln(saved.checkpoint_signature, checkpoint_signature) && ...
        height(saved.detailed_table) == n_rows;

    if compatible
        detailed_table = saved.detailed_table;
        candidate_complete = saved.candidate_complete;
        fprintf("Resuming: %d of %d candidates already complete.\n\n", ...
            nnz(candidate_complete), n_candidates);
    else
        warning("Existing checkpoint is incompatible and will be ignored.");
    end
end

%% Run the six combinations with identical conditions

overall_clock = tic;

for i_candidate = 1:n_candidates
    if candidate_complete(i_candidate)
        continue
    end

    candidate_clock = tic;
    tau_s = candidate_tau_s(i_candidate);
    delay_s = candidate_delay_s(i_candidate);
    drag_xy = candidate_drag_xy(i_candidate);

    fprintf("[%d/%d] %s | tau %.1f ms | delay %.1f ms | d_q,xy %.4f\n", ...
        i_candidate, n_candidates, candidate_labels(i_candidate), ...
        1000*tau_s, 1000*delay_s, drag_xy);

    P = P_base;
    P.actuators.servo_time_constant = tau_s;
    P.feedback.state_delay = delay_s;
    P.aerodynamics.quadratic_drag_B(1,1) = drag_xy;
    P.aerodynamics.quadratic_drag_B(2,2) = drag_xy;

    for i_run = 1:n_runs
        row = (i_candidate - 1)*n_runs + i_run;
        test = tests(runs(i_run).test_index);
        configuration = runs(i_run).configuration;
        target = test.(char(configuration));

        scenario = struct();
        scenario.type = test.trajectory;
        scenario.vmax = test.vmax;
        scenario.amax = test.amax;
        scenario.cycles = P.trajectory.default_cycles;
        scenario.start_from_rest = true;

        fprintf("    %-11s %-6s", test.id, configuration);
        run_clock = tic;
        result = simulate_COMET(P, configuration, scenario);
        metrics = result.metrics;

        detailed_table.Candidate(row) = i_candidate;
        detailed_table.Candidate_label(row) = candidate_labels(i_candidate);
        detailed_table.Servo_tau_s(row) = tau_s;
        [applied_delay_s, applied_delay_samples] = ...
            read_applied_feedback_delay(result, delay_s, controller_dt, ...
                simulator_path);
        detailed_table.Feedback_delay_s(row) = applied_delay_s;
        detailed_table.Feedback_delay_samples(row) = applied_delay_samples;
        detailed_table.Quadratic_drag_xy(row) = drag_xy;
        detailed_table.Case(row) = test.id;
        detailed_table.Trajectory(row) = test.trajectory;
        detailed_table.Configuration(row) = configuration;
        detailed_table.Vmax_mps(row) = test.vmax;
        detailed_table.Experimental_available(row) = target.available;
        detailed_table.Experimental_completed(row) = target.completed;
        detailed_table.Expected_failure(row) = ...
            target.available && ~target.completed;

        detailed_table.RMSE_sim_m(row) = metrics.rmse_position;
        detailed_table.Max_error_sim_m(row) = metrics.max_position_error;
        detailed_table.Velocity_RMSE_sim_mps(row) = metrics.rmse_velocity;
        detailed_table.Roll_peak_sim(row) = ...
            max(abs(metrics.roll_command_range));
        detailed_table.Pitch_peak_sim(row) = ...
            max(abs(metrics.pitch_command_range));
        detailed_table.Cyclic_position_saturation_percent(row) = ...
            100*metrics.cyclic_position_saturation_fraction;
        detailed_table.Rate_saturation_percent(row) = ...
            100*metrics.actuator_rate_saturation_fraction;
        detailed_table.Any_saturation_percent(row) = ...
            100*metrics.any_saturation_fraction;
        detailed_table.Peak_cyclic_actual_percent(row) = ...
            peak_cyclic_use(result);
        detailed_table.Lost_control(row) = metrics.lost_control;

        if target.available && target.completed
            detailed_table.RMSE_exp_m(row) = target.rmse_position;
            detailed_table.Max_error_exp_m(row) = ...
                target.max_position_error;
            detailed_table.Roll_peak_exp(row) = ...
                max(abs(target.roll_command_range));
            detailed_table.Pitch_peak_exp(row) = ...
                max(abs(target.pitch_command_range));

            detailed_table.RMSE_delta_percent(row) = percent_mismatch( ...
                detailed_table.RMSE_sim_m(row), ...
                detailed_table.RMSE_exp_m(row));
            detailed_table.Max_error_delta_percent(row) = percent_mismatch( ...
                detailed_table.Max_error_sim_m(row), ...
                detailed_table.Max_error_exp_m(row));
            detailed_table.Roll_peak_delta_percent(row) = percent_mismatch( ...
                detailed_table.Roll_peak_sim(row), ...
                detailed_table.Roll_peak_exp(row));
            detailed_table.Pitch_peak_delta_percent(row) = percent_mismatch( ...
                detailed_table.Pitch_peak_sim(row), ...
                detailed_table.Pitch_peak_exp(row));

            detailed_table.Both_position_metrics_within_10_percent(row) = ...
                abs(detailed_table.RMSE_delta_percent(row)) <= ...
                    acceptable_mismatch_percent && ...
                abs(detailed_table.Max_error_delta_percent(row)) <= ...
                    acceptable_mismatch_percent && ...
                ~detailed_table.Lost_control(row);
        end

        fprintf(" | RMSE %.3f m | max %.3f m | sat %.1f %% | %.1f s\n", ...
            metrics.rmse_position, metrics.max_position_error, ...
            max(100*metrics.cyclic_position_saturation_fraction, ...
                100*metrics.actuator_rate_saturation_fraction), ...
            toc(run_clock));
        drawnow
        clear result metrics
    end

    candidate_complete(i_candidate) = true;

    if save_results
        save(checkpoint_path, "detailed_table", "candidate_complete", ...
            "candidate_tau_s", "candidate_delay_s", ...
            "candidate_drag_xy", "candidate_labels", ...
            "checkpoint_signature", "acceptable_mismatch_percent", ...
            "completed_case_saturation_guardrail_percent");
    end

    completed_count = nnz(candidate_complete);
    mean_candidate_time = toc(overall_clock)/max(completed_count,1);
    remaining_minutes = ...
        mean_candidate_time*(n_candidates - completed_count)/60;
    fprintf("    completed in %.1f s | estimated remaining %.1f min\n\n", ...
        toc(candidate_clock), remaining_minutes);
end

%% Position-only score and ranking

candidate_summary = summarize_candidates(detailed_table, candidates, ...
    acceptable_mismatch_percent, ...
    completed_case_saturation_guardrail_percent);

selection_matrix = [ ...
    candidate_summary.Unexpected_lost_control_count, ...
    candidate_summary.Completed_saturation_guardrail_violations, ...
    candidate_summary.Position_score_percent];
[~, rank_order] = sortrows(selection_matrix, [1, 2, 3]);
ranked_summary = candidate_summary(rank_order,:);
ranked_summary.Rank = (1:n_candidates).';
ranked_summary = movevars(ranked_summary, "Rank", "Before", "Candidate");

disp("POSITION-METRIC CANDIDATE RANKING")
disp(ranked_summary)

best = ranked_summary(1,:);
best_is_admissible = best.Unexpected_lost_control_count == 0 && ...
    best.Completed_saturation_guardrail_violations == 0;

if best_is_admissible
    selection_description = "best admissible candidate";
else
    selection_description = ...
        "least-bad diagnostic candidate (no candidate is fully admissible)";
    warning(["No candidate satisfies every completed-case guardrail. " + ...
        "The first candidate must not yet be treated as calibrated."]);
end

fprintf("\n%s: %s\n", upper(selection_description), best.Label);
fprintf("tau %.1f ms | delay %.1f ms | d_q,xy %.4f\n", ...
    best.Servo_tau_ms, best.Feedback_delay_ms, best.Quadratic_drag_xy);
fprintf("Position score %.1f %% | worst mismatch %.1f %% | within 10%% %.1f %%\n", ...
    best.Position_score_percent, ...
    best.Worst_position_mismatch_percent, ...
    best.Position_metrics_within_10_percent);

fprintf("\nAssignments to test this candidate in RUN_COMET/comet_parameters:\n");
fprintf("P.actuators.servo_time_constant = %.6f;\n", ...
    best.Servo_tau_ms/1000);
fprintf("P.feedback.state_delay = %.6f;\n", ...
    best.Feedback_delay_ms/1000);
fprintf("P.aerodynamics.quadratic_drag_B(1,1) = %.6f;\n", ...
    best.Quadratic_drag_xy);
fprintf("P.aerodynamics.quadratic_drag_B(2,2) = %.6f;\n\n", ...
    best.Quadratic_drag_xy);

%% Compact plots

compact_position_path = fullfile(results_folder, ...
    "COMET_targeted_compact_position_metrics.png");
plot_compact_position_metrics(detailed_table, candidates, ...
    acceptable_mismatch_percent, ...
    completed_case_saturation_guardrail_percent, ...
    compact_position_path, save_results);

ranking_path = fullfile(results_folder, ...
    "COMET_targeted_candidate_ranking.png");
plot_candidate_ranking(candidate_summary, best.Candidate, ...
    acceptable_mismatch_percent, ranking_path, save_results);

%% Save numerical outputs

selected_parameters = struct();
selected_parameters.label = best.Label;
selected_parameters.servo_time_constant_s = best.Servo_tau_ms/1000;
selected_parameters.feedback_delay_s = best.Feedback_delay_ms/1000;
selected_parameters.quadratic_drag_xy = best.Quadratic_drag_xy;
selected_parameters.position_score_percent = best.Position_score_percent;
selected_parameters.admissible = best_is_admissible;

if save_results
    writetable(candidates, fullfile(results_folder, ...
        "COMET_targeted_candidate_definitions.csv"));
    writetable(detailed_table, fullfile(results_folder, ...
        "COMET_targeted_detailed_results.csv"));
    writetable(candidate_summary, fullfile(results_folder, ...
        "COMET_targeted_candidate_summary.csv"));
    writetable(ranked_summary, fullfile(results_folder, ...
        "COMET_targeted_ranked_summary.csv"));

    save(fullfile(results_folder, "COMET_targeted_results.mat"), ...
        "detailed_table", "candidate_summary", "ranked_summary", ...
        "candidates", "selected_parameters", ...
        "acceptable_mismatch_percent", ...
        "completed_case_saturation_guardrail_percent");
    save(fullfile(results_folder, "COMET_selected_parameters.mat"), ...
        "selected_parameters");

    fprintf("Results saved in:\n%s\n", results_folder);
end


%% Local functions

function runs = make_run_list(tests)

runs = repmat(struct("test_index", 0, "configuration", ""), 0, 1);
for i_test = 1:numel(tests)
    for configuration = ["single", "dual"]
        entry.test_index = i_test;
        entry.configuration = configuration;
        runs(end+1,1) = entry; %#ok<AGROW>
    end
end

end


function [applied_delay_s, applied_delay_samples] = ...
    read_applied_feedback_delay(result, requested_delay_s, ...
        controller_dt, simulator_path)

if isfield(result, "notes") && ...
        isfield(result.notes, "feedback_delay_applied_s") && ...
        isfield(result.notes, "feedback_delay_samples")
    applied_delay_s = result.notes.feedback_delay_applied_s;
    applied_delay_samples = result.notes.feedback_delay_samples;
elseif isfield(result, "parameters") && ...
        isfield(result.parameters, "feedback") && ...
        isfield(result.parameters.feedback, "state_delay_applied") && ...
        isfield(result.parameters.feedback, "state_delay_samples")
    applied_delay_s = result.parameters.feedback.state_delay_applied;
    applied_delay_samples = result.parameters.feedback.state_delay_samples;
else
    error("CALIBRATE_COMET_MULTIPARAMETER:DelayModelMissing", ...
        ["The active simulate_COMET.m does not expose feedback-delay " + ...
         "diagnostics. MATLAB is using: %s"], simulator_path);
end

expected_samples = round(requested_delay_s/controller_dt);
expected_delay_s = expected_samples*controller_dt;

if applied_delay_samples ~= expected_samples || ...
        abs(applied_delay_s - expected_delay_s) > ...
        10*eps(max(1,expected_delay_s))
    error("CALIBRATE_COMET_MULTIPARAMETER:DelayNotApplied", ...
        ["Requested %.1f ms (%d samples), obtained %.1f ms (%d samples). " + ...
         "Active file: %s"], ...
        1000*expected_delay_s, expected_samples, ...
        1000*applied_delay_s, applied_delay_samples, simulator_path);
end

end


function signature = make_checkpoint_signature(P, script_folder)

parameter_file = dir(fullfile(script_folder, "comet_parameters.m"));
simulation_file = dir(fullfile(script_folder, "simulate_COMET.m"));

signature.parameter_version = P.meta.version;
signature.K_cyc = P.aerodynamics.cyclic_effectiveness_scale;
signature.controller_dt = P.simulation.controller_sample_time;
signature.parameter_bytes = parameter_file.bytes;
signature.parameter_timestamp = parameter_file.datenum;
signature.simulation_bytes = simulation_file.bytes;
signature.simulation_timestamp = simulation_file.datenum;

end


function T = empty_detailed_table(n)

Candidate = nan(n,1);
Candidate_label = strings(n,1);
Servo_tau_s = nan(n,1);
Feedback_delay_s = nan(n,1);
Feedback_delay_samples = nan(n,1);
Quadratic_drag_xy = nan(n,1);
Case = strings(n,1);
Trajectory = strings(n,1);
Configuration = strings(n,1);
Vmax_mps = nan(n,1);
Experimental_available = false(n,1);
Experimental_completed = false(n,1);
Expected_failure = false(n,1);
RMSE_sim_m = nan(n,1);
RMSE_exp_m = nan(n,1);
RMSE_delta_percent = nan(n,1);
Max_error_sim_m = nan(n,1);
Max_error_exp_m = nan(n,1);
Max_error_delta_percent = nan(n,1);
Velocity_RMSE_sim_mps = nan(n,1);
Roll_peak_sim = nan(n,1);
Roll_peak_exp = nan(n,1);
Roll_peak_delta_percent = nan(n,1);
Pitch_peak_sim = nan(n,1);
Pitch_peak_exp = nan(n,1);
Pitch_peak_delta_percent = nan(n,1);
Cyclic_position_saturation_percent = nan(n,1);
Rate_saturation_percent = nan(n,1);
Any_saturation_percent = nan(n,1);
Peak_cyclic_actual_percent = nan(n,1);
Lost_control = false(n,1);
Both_position_metrics_within_10_percent = false(n,1);

T = table(Candidate, Candidate_label, Servo_tau_s, Feedback_delay_s, ...
    Feedback_delay_samples, Quadratic_drag_xy, Case, Trajectory, ...
    Configuration, Vmax_mps, Experimental_available, ...
    Experimental_completed, Expected_failure, RMSE_sim_m, RMSE_exp_m, ...
    RMSE_delta_percent, Max_error_sim_m, Max_error_exp_m, ...
    Max_error_delta_percent, Velocity_RMSE_sim_mps, ...
    Roll_peak_sim, Roll_peak_exp, Roll_peak_delta_percent, ...
    Pitch_peak_sim, Pitch_peak_exp, Pitch_peak_delta_percent, ...
    Cyclic_position_saturation_percent, Rate_saturation_percent, ...
    Any_saturation_percent, Peak_cyclic_actual_percent, ...
    Lost_control, Both_position_metrics_within_10_percent);

end


function value = percent_mismatch(simulated, experimental)

if ~isfinite(experimental) || abs(experimental) < eps
    value = NaN;
else
    value = 100*(simulated - experimental)/abs(experimental);
end

end


function value = mean_finite(values)

values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = mean(values);
end

end


function value = max_finite(values)

values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = max(values);
end

end


function S = summarize_candidates(T, candidates, tolerance, guardrail)

n_candidates = height(candidates);
Candidate = candidates.Candidate;
Label = candidates.Label;
Servo_tau_ms = candidates.Servo_tau_ms;
Feedback_delay_ms = candidates.Feedback_delay_ms;
Quadratic_drag_xy = candidates.Quadratic_drag_xy;
RMSE_MAPE_percent = nan(n_candidates,1);
Max_error_MAPE_percent = nan(n_candidates,1);
Position_score_percent = nan(n_candidates,1);
Worst_position_mismatch_percent = nan(n_candidates,1);
Position_metrics_within_10_percent = nan(n_candidates,1);
Unexpected_lost_control_count = zeros(n_candidates,1);
Completed_saturation_guardrail_violations = zeros(n_candidates,1);
Expected_failure_missed_count = zeros(n_candidates,1);

for i_candidate = 1:n_candidates
    rows = T.Candidate == i_candidate;
    completed = rows & T.Experimental_completed;
    expected_failure = rows & T.Expected_failure;

    rmse_mismatch = abs(T.RMSE_delta_percent(completed));
    max_mismatch = abs(T.Max_error_delta_percent(completed));
    all_position_mismatches = [rmse_mismatch; max_mismatch];

    RMSE_MAPE_percent(i_candidate) = mean_finite(rmse_mismatch);
    Max_error_MAPE_percent(i_candidate) = mean_finite(max_mismatch);
    Position_score_percent(i_candidate) = ...
        mean_finite(all_position_mismatches);
    Worst_position_mismatch_percent(i_candidate) = ...
        max_finite(all_position_mismatches);
    Position_metrics_within_10_percent(i_candidate) = ...
        100*mean_finite(all_position_mismatches <= tolerance);

    Unexpected_lost_control_count(i_candidate) = ...
        sum(T.Lost_control(completed));

    completed_saturation = max( ...
        T.Cyclic_position_saturation_percent(completed), ...
        T.Rate_saturation_percent(completed));
    Completed_saturation_guardrail_violations(i_candidate) = ...
        sum(completed_saturation > guardrail);

    reproduced_failure = T.Lost_control(expected_failure) & ...
        (T.Cyclic_position_saturation_percent(expected_failure) > 0 | ...
         T.Rate_saturation_percent(expected_failure) > 0);
    Expected_failure_missed_count(i_candidate) = sum(~reproduced_failure);
end

S = table(Candidate, Label, Servo_tau_ms, Feedback_delay_ms, ...
    Quadratic_drag_xy, RMSE_MAPE_percent, Max_error_MAPE_percent, ...
    Position_score_percent, Worst_position_mismatch_percent, ...
    Position_metrics_within_10_percent, ...
    Unexpected_lost_control_count, ...
    Completed_saturation_guardrail_violations, ...
    Expected_failure_missed_count);

end


function plot_compact_position_metrics(T, candidates, tolerance, ...
    guardrail, output_path, save_results)

case_order = unique(T.Case, "stable");
n_cases = numel(case_order);
n_candidates = height(candidates);
x = (1:n_candidates).';

single_sim_color = [0.000, 0.447, 0.741];
single_exp_color = [0.301, 0.745, 0.933];
dual_sim_color = [0.850, 0.325, 0.098];
dual_exp_color = [0.929, 0.694, 0.125];

metric_sim = ["RMSE_sim_m", "Max_error_sim_m"];
metric_exp = ["RMSE_exp_m", "Max_error_exp_m"];
metric_ylabel = ["Position RMSE [m]", "Maximum position error [m]"];

figure("Color", "w", "Name", "COMET targeted position comparison", ...
    "Position", [20, 40, 2350, 1050]);
tiledlayout(2, n_cases, "TileSpacing", "compact", "Padding", "compact");

for i_metric = 1:2
    for i_case = 1:n_cases
        ax = nexttile((i_metric - 1)*n_cases + i_case);
        hold(ax, "on")

        case_name = case_order(i_case);
        single = T(T.Case == case_name & ...
            T.Configuration == "single", :);
        dual = T(T.Case == case_name & ...
            T.Configuration == "dual", :);
        single = sortrows(single, "Candidate");
        dual = sortrows(dual, "Candidate");

        y_single = single.(char(metric_sim(i_metric)));
        y_dual = dual.(char(metric_sim(i_metric)));
        exp_single = first_finite(single.(char(metric_exp(i_metric))));
        exp_dual = first_finite(dual.(char(metric_exp(i_metric))));

        if isfinite(exp_single)
            lower = exp_single*(1 - tolerance/100);
            upper = exp_single*(1 + tolerance/100);
            fill(ax, [0.5, n_candidates+0.5, n_candidates+0.5, 0.5], ...
                [lower, lower, upper, upper], single_exp_color, ...
                "FaceAlpha", 0.10, "EdgeColor", "none", ...
                "HandleVisibility", "off");
            h_single_exp = plot(ax, [0.5, n_candidates+0.5], ...
                [exp_single, exp_single], "--", ...
                "Color", single_exp_color, "LineWidth", 1.8);
        else
            h_single_exp = plot(ax, nan, nan, "--", ...
                "Color", single_exp_color, "LineWidth", 1.8);
        end

        if isfinite(exp_dual)
            lower = exp_dual*(1 - tolerance/100);
            upper = exp_dual*(1 + tolerance/100);
            fill(ax, [0.5, n_candidates+0.5, n_candidates+0.5, 0.5], ...
                [lower, lower, upper, upper], dual_exp_color, ...
                "FaceAlpha", 0.10, "EdgeColor", "none", ...
                "HandleVisibility", "off");
            h_dual_exp = plot(ax, [0.5, n_candidates+0.5], ...
                [exp_dual, exp_dual], "--", ...
                "Color", dual_exp_color, "LineWidth", 1.8);
        else
            h_dual_exp = plot(ax, nan, nan, "--", ...
                "Color", dual_exp_color, "LineWidth", 1.8);
        end

        h_single_sim = plot(ax, x, y_single, "o-", ...
            "Color", single_sim_color, "MarkerFaceColor", "w", ...
            "LineWidth", 1.8, "MarkerSize", 6);
        h_dual_sim = plot(ax, x, y_dual, "s-", ...
            "Color", dual_sim_color, "MarkerFaceColor", "w", ...
            "LineWidth", 1.8, "MarkerSize", 6);

        single_bad_sat = max( ...
            single.Cyclic_position_saturation_percent, ...
            single.Rate_saturation_percent) > guardrail;
        dual_bad_sat = max( ...
            dual.Cyclic_position_saturation_percent, ...
            dual.Rate_saturation_percent) > guardrail;
        bad_x = [x(single_bad_sat); x(dual_bad_sat)];
        bad_y = [y_single(single_bad_sat); y_dual(dual_bad_sat)];
        h_saturation = plot(ax, bad_x, bad_y, "kx", ...
            "LineWidth", 1.5, "MarkerSize", 9);

        lost_x = [x(single.Lost_control); x(dual.Lost_control)];
        lost_y = [y_single(single.Lost_control); y_dual(dual.Lost_control)];
        h_lost = plot(ax, lost_x, lost_y, "rp", ...
            "MarkerFaceColor", "r", "MarkerSize", 10);

        if i_metric == 1
            title(ax, strrep(case_name, "_", " "))
        end
        ylabel(ax, metric_ylabel(i_metric))
        xticks(ax, x)
        xticklabels(ax, candidates.Label)
        xlim(ax, [0.5, n_candidates+0.5])
        grid(ax, "on")

        finite_values = [y_single(isfinite(y_single)); ...
            y_dual(isfinite(y_dual)); exp_single; exp_dual];
        finite_values = finite_values(isfinite(finite_values));
        if isempty(finite_values)
            y_upper = 1;
        else
            y_upper = max(finite_values);
        end
        ylim(ax, [0, max(1.08*y_upper, eps)])

        if ~single.Experimental_completed(1)
            if single.Expected_failure(1)
                note = "Single experiment: loss of control";
            else
                note = "Single experiment: not reported";
            end
            text(ax, 0.03, 0.96, note, "Units", "normalized", ...
                "VerticalAlignment", "top", "Color", [0.55, 0.10, 0.10], ...
                "FontSize", 7.5, "FontWeight", "bold");
        end

        if i_metric == 2
            xlabel(ax, "Candidate combination")
        end

        if i_metric == 1 && i_case == 1
            legend(ax, [h_single_sim, h_single_exp, ...
                h_dual_sim, h_dual_exp, h_saturation, h_lost], ...
                {"Simulation single", "Experiment single", ...
                 "Simulation dual", "Experiment dual", ...
                 "Saturation > 5%", "Lost control"}, ...
                "Location", "best", "FontSize", 7);
        end
    end
end

mapping = strings(n_candidates,1);
for i = 1:n_candidates
    mapping(i) = sprintf("%s: delay %.0f ms, d_{q,xy}=%.4f", ...
        candidates.Label(i), candidates.Feedback_delay_ms(i), ...
        candidates.Quadratic_drag_xy(i));
end

sgtitle(sprintf(["Targeted delay/drag comparison; tau = %.0f ms; " + ...
    "shaded regions = experiment +/-%.0f%%"], ...
    candidates.Servo_tau_ms(1), tolerance));
annotation(gcf, "textbox", [0.04, 0.002, 0.92, 0.035], ...
    "String", strjoin(mapping, "    "), "EdgeColor", "none", ...
    "HorizontalAlignment", "center", "FontSize", 8, ...
    "Interpreter", "tex");

if save_results
    exportgraphics(gcf, output_path, "Resolution", 180);
end

end


function value = first_finite(values)

values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = values(1);
end

end


function plot_candidate_ranking(S, selected_candidate, tolerance, ...
    output_path, save_results)

x = S.Candidate;

figure("Color", "w", "Name", "COMET targeted candidate ranking", ...
    "Position", [100, 100, 1500, 650]);
tiledlayout(1,2, "TileSpacing", "compact", "Padding", "compact");

nexttile
b = bar(x, [S.RMSE_MAPE_percent, S.Max_error_MAPE_percent], "grouped");
b(1).FaceColor = [0.000, 0.447, 0.741];
b(2).FaceColor = [0.850, 0.325, 0.098];
hold on
yline(tolerance, "k--", sprintf("%.0f%% target", tolerance), ...
    "LineWidth", 1.3);
selected_height = max([S.RMSE_MAPE_percent(selected_candidate), ...
    S.Max_error_MAPE_percent(selected_candidate)]);
plot(selected_candidate, selected_height*1.08, ...
    "kp", "MarkerFaceColor", [0.929, 0.694, 0.125], ...
    "MarkerSize", 13);
xticks(x)
xticklabels(S.Label)
xlabel("Candidate combination")
ylabel("Mean absolute mismatch [%]")
title("Position-metric mismatch")
legend("RMSE", "Maximum error", "Location", "best")
grid on

nexttile
guardrail_values = [S.Unexpected_lost_control_count, ...
    S.Completed_saturation_guardrail_violations, ...
    S.Expected_failure_missed_count];
b = bar(x, guardrail_values, "grouped");
b(1).FaceColor = [0.635, 0.078, 0.184];
b(2).FaceColor = [0.250, 0.250, 0.250];
b(3).FaceColor = [0.466, 0.674, 0.188];
xticks(x)
xticklabels(S.Label)
xlabel("Candidate combination")
ylabel("Number of cases")
title("Physical guardrails")
legend("Unexpected loss", "Completed-case saturation >5%", ...
    "Expected failure missed", "Location", "best")
grid on

sgtitle("Targeted candidate ranking; star = selected candidate")

if save_results
    exportgraphics(gcf, output_path, "Resolution", 180);
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
