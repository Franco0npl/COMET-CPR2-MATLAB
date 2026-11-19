function out = simulate_COMET(P, configuration, scenario)
%SIMULATE_COMET Simulate the COMET vehicle with a fixed-step RK4 scheme.
%
%   OUT = SIMULATE_COMET(P, CONFIGURATION, SCENARIO) runs one nonlinear
%   trajectory-tracking simulation. The same controller, gains and mixer are
%   used for both configurations. The only architectural change is the cyclic
%   mask:
%
%       configuration = "single"  -> upper cyclic disabled, lower active
%       configuration = "dual"    -> upper and lower cyclic active
%
%   P is returned by comet_parameters(). SCENARIO may be a string
%   ("hover", "circle" or "eight") or a structure. Typical use is:
%
%       P = comet_parameters();
%       S.type   = "circle";
%       S.vmax   = 4.0;
%       S.cycles = 2;
%
%       single = simulate_COMET(P, "single", S);
%       dual   = simulate_COMET(P, "dual",   S);
%
%   Main state convention (scalar-first quaternion):
%
%       x = [p_I; v_I; q_BI; omega_B; actuator_state]
%
%   where q_BI and R(q_BI) rotate body-frame vectors into the inertial frame.
%   actuator_state is ordered as
%
%       [T_up; T_dw; ail_up; ele_up; ail_dw; ele_dw].
%
%   The controller is evaluated once every
%   P.simulation.controller_sample_time seconds. Its output is held constant
%   over that complete control interval. RK4 can divide the interval into
%   P.simulation.integration_substeps without changing the controller rate.
%   Actuator states remain continuous and are integrated inside every RK4
%   substep. Numerical error derivatives are low-pass filtered before
%   multiplication by Kd to avoid sample-to-sample cyclic chattering.
%
%   This is a reduced-order COMET model. In particular, the fields named
%   cyclic_arm_* in P are interpreted here as residual hub/flapping-moment
%   coefficients. The geometric moment r_hub x F_rotor is computed separately
%   and is therefore not hidden inside those coefficients.

if nargin < 1 || isempty(P)
    P = comet_parameters();
end
if nargin < 2 || isempty(configuration)
    configuration = "dual";
end
if nargin < 3 || isempty(scenario)
    scenario = struct();
end

I = state_indices();
C = select_configuration(P, configuration);
S = complete_scenario(P, scenario);

controller_dt = P.simulation.controller_sample_time;
integration_substeps = P.simulation.integration_substeps;
integration_dt = controller_dt/integration_substeps;

% The delay acts only on the state seen by the feedback controller. Delays
% that are not integer multiples of the controller sample time are obtained
% by linear interpolation between the two surrounding stored states. This
% preserves the requested delay instead of silently rounding it.
feedback_delay_requested = 0.0;
if isfield(P, "feedback") && isfield(P.feedback, "state_delay")
    feedback_delay_requested = P.feedback.state_delay;
end
if ~isscalar(feedback_delay_requested) || ...
        ~isfinite(feedback_delay_requested) || ...
        feedback_delay_requested < 0
    error('simulate_COMET:InvalidFeedbackDelay', ...
        "P.feedback.state_delay must be a finite nonnegative scalar.");
end
feedback_delay_samples = feedback_delay_requested/controller_dt;
feedback_delay_floor_samples = floor(feedback_delay_samples);
feedback_delay_fraction = ...
    feedback_delay_samples - feedback_delay_floor_samples;

% Remove floating-point residue for delays that are exact sample multiples.
delay_tolerance = 1e-12;
if feedback_delay_fraction < delay_tolerance
    feedback_delay_fraction = 0.0;
elseif 1.0 - feedback_delay_fraction < delay_tolerance
    feedback_delay_floor_samples = feedback_delay_floor_samples + 1;
    feedback_delay_fraction = 0.0;
end

feedback_delay_applied = feedback_delay_requested;
P.feedback.state_delay_requested = feedback_delay_requested;
P.feedback.state_delay_applied = feedback_delay_applied;
P.feedback.state_delay_samples = feedback_delay_samples;
P.feedback.state_delay_floor_samples = feedback_delay_floor_samples;
P.feedback.state_delay_fraction = feedback_delay_fraction;

n_intervals = ceil(S.final_time / controller_dt);
t = (0:n_intervals).' * controller_dt;
N = numel(t);
S.simulated_final_time = t(end);

x = initial_state(P, S, I);
memory = initial_controller_memory();
log = preallocate_log(N, I.n_actuators);
feedback_state_history = nan(I.n_states, N);

last_index = N;
termination_reason = "completed_time_horizon";
numerical_divergence = false;
stopped_for_lost_control = false;

for k = 1:N
    tk = t(k);
    ref = trajectory_reference(P, S, tk);

    feedback_state_history(:,k) = x;
    feedback_recent_index = max(1, k - feedback_delay_floor_samples);
    feedback_recent = feedback_state_history(:,feedback_recent_index);

    if feedback_delay_fraction == 0.0
        feedback_state = feedback_recent;
    else
        feedback_older_index = max(1, feedback_recent_index - 1);
        feedback_older = feedback_state_history(:,feedback_older_index);
        feedback_state = (1.0 - feedback_delay_fraction)*feedback_recent + ...
            feedback_delay_fraction*feedback_older;

        % Quaternions q and -q represent the same attitude. Align their
        % signs before interpolation, then restore unit norm.
        q_recent = feedback_recent(I.quaternion);
        q_older = feedback_older(I.quaternion);
        if dot(q_recent, q_older) < 0.0
            q_older = -q_older;
        end
        q_feedback = (1.0 - feedback_delay_fraction)*q_recent + ...
            feedback_delay_fraction*q_older;
        feedback_state(I.quaternion) = normalize_quaternion(q_feedback);
    end

    [control, memory_sample] = controller_sample( ...
        P, feedback_state, ref, memory, controller_dt, I);
    allocation = allocate_actuators(P, C, control);
    rate_status = actuator_rate(P, C, x(I.actuators), ...
        allocation.command_limited);

    [rotor_wrench, rotor_detail] = rotor_wrench_body(P, C, ...
        x(I.actuators));
    [environment, environment_detail] = environmental_wrench(P, tk, x, I);

    log = store_sample(log, k, x, ref, control, allocation, ...
        rate_status, rotor_wrench, rotor_detail, environment, ...
        environment_detail, P, I);

    position_error_norm = norm(ref.position_I - x(I.position));
    lost_control_now = tk >= S.metric_start_time && ...
        position_error_norm > P.simulation.lost_control_position_error;

    if lost_control_now && P.simulation.stop_on_lost_control
        last_index = k;
        termination_reason = "lost_control_threshold";
        stopped_for_lost_control = true;
        break
    end

    if k == N
        break
    end

    % The actuator command is fixed over the complete controller interval.
    % Only the numerical plant integration is refined by the substep loop.
    dynamics = @(tau, state) state_derivative( ...
        tau, state, allocation.command_limited, P, C, I);
    x_candidate = x;
    projection_ok = true;
    state_ok = true;
    state_failure = "state_valid";
    failure_time = tk + controller_dt;

    for substep = 1:integration_substeps
        substep_time = tk + (substep - 1)*integration_dt;
        x_candidate = rk4_step( ...
            dynamics, substep_time, x_candidate, integration_dt);
        [x_candidate, projection_ok] = ...
            project_state(P, C, x_candidate, I);
        [state_ok, state_failure] = ...
            validate_dynamic_state(P, x_candidate, I);
        failure_time = substep_time + integration_dt;
        if ~(projection_ok && state_ok)
            break
        end
    end

    if ~(projection_ok && state_ok)
        last_index = k;
        termination_reason = state_failure;
        numerical_divergence = true;
        if ~P.simulation.stop_on_numerical_divergence
            error('simulate_COMET:NumericalDivergence', ...
                "Numerical divergence at t = %.6f s: %s", ...
                failure_time, char(state_failure));
        end
        break
    end

    saturation_active = any(allocation.position_saturated) || ...
        any(rate_status.saturated);
    memory = update_controller_memory(P, memory_sample, control, ...
        saturation_active, controller_dt);
    x = x_candidate;
end

out = assemble_output(P, C, S, t, log, last_index, ...
    termination_reason, numerical_divergence, stopped_for_lost_control);

end


%% State layout

function I = state_indices()

I.position = 1:3;
I.velocity = 4:6;
I.quaternion = 7:10;
I.omega = 11:13;
I.actuators = 14:19;
I.n_states = 19;
I.n_actuators = 6;

end


%% Configuration and scenario preparation

function C = select_configuration(P, configuration)

name = lower(strtrim(string(configuration)));

switch name
    case {"single", "single-lower", "single_lower"}
        C.name = P.configuration.single.name;
        C.cyclic_mask = P.configuration.single.cyclic_mask(:);
        C.is_dual = false;

    case "dual"
        C.name = P.configuration.dual.name;
        C.cyclic_mask = P.configuration.dual.cyclic_mask(:);
        C.is_dual = true;

    otherwise
        error('simulate_COMET:UnknownConfiguration', ...
            "Unknown configuration '%s'. Use 'single' or 'dual'.", ...
            char(name));
end

assert(numel(C.cyclic_mask) == 4, ...
    "The cyclic mask must contain four entries.");
assert(all(C.cyclic_mask == 0 | C.cyclic_mask == 1), ...
    "The cyclic mask must contain only zeros and ones.");

end


function S = complete_scenario(P, scenario)

if ischar(scenario) || (isstring(scenario) && isscalar(scenario))
    input_scenario.type = string(scenario);
elseif isstruct(scenario)
    input_scenario = scenario;
else
    error('simulate_COMET:InvalidScenario', ...
        "SCENARIO must be a trajectory name or a structure.");
end

if ~isfield(input_scenario, "type")
    if isfield(input_scenario, "trajectory")
        input_scenario.type = input_scenario.trajectory;
    else
        input_scenario.type = "circle";
    end
end

S.type = lower(strtrim(string(input_scenario.type)));
if S.type == "figure_eight" || S.type == "figure-eight"
    S.type = "eight";
end

S.height = get_option(input_scenario, "height", ...
    P.trajectory.hover_height);
S.yaw = get_option(input_scenario, "yaw", ...
    P.trajectory.default_yaw);
S.cycles = get_option(input_scenario, "cycles", ...
    P.trajectory.default_cycles);
S.start_from_rest = logical(get_option(input_scenario, ...
    "start_from_rest", P.trajectory.start_from_rest));

switch S.type
    case "hover"
        S.vmax = 0;
        S.amax = 0;
        S.phase_rate = 0;
        S.target_phase = 0;
        S.entry_time = 0;
        S.final_time = get_option(input_scenario, "duration", 10.0);

    case "circle"
        S.vmax = get_option(input_scenario, "vmax", 3.0);
        S.radius = get_option(input_scenario, "radius", ...
            P.trajectory.circle_radius);
        require_positive(S.vmax, "scenario.vmax");
        require_positive(S.radius, "scenario.radius");
        require_positive(S.cycles, "scenario.cycles");

        S.phase_rate = S.vmax / S.radius;
        S.target_phase = 2*pi*S.cycles;
        theoretical_amax = S.vmax^2 / S.radius;
        S.amax = get_option(input_scenario, "amax", theoretical_amax);
        S.entry_time = choose_entry_time(P, input_scenario, S);
        S.final_time = choose_final_time(input_scenario, S);

    case "eight"
        S.vmax = get_option(input_scenario, "vmax", 3.5);
        S.A = get_option(input_scenario, "A", P.trajectory.eight_A);
        S.B = get_option(input_scenario, "B", P.trajectory.eight_B);
        require_positive(S.vmax, "scenario.vmax");
        require_positive(S.A, "scenario.A");
        require_positive(S.B, "scenario.B");
        require_positive(S.cycles, "scenario.cycles");

        % For x=A sin(theta), y=B sin(2 theta), the maximum value of
        % ||dr/dtheta|| occurs at theta=0 and is sqrt(A^2+(2B)^2).
        speed_scale = sqrt(S.A^2 + (2*S.B)^2);
        S.phase_rate = S.vmax / speed_scale;
        S.target_phase = 2*pi*S.cycles;

        phase_grid = linspace(0, 2*pi, 4001);
        ax = -S.A * sin(phase_grid) * S.phase_rate^2;
        ay = -4*S.B * sin(2*phase_grid) * S.phase_rate^2;
        theoretical_amax = max(hypot(ax, ay));
        S.amax = get_option(input_scenario, "amax", theoretical_amax);
        S.entry_time = choose_entry_time(P, input_scenario, S);
        S.final_time = choose_final_time(input_scenario, S);

    otherwise
        error('simulate_COMET:UnknownTrajectory', ...
            "Unknown trajectory '%s'. Use hover, circle or eight.", ...
            char(S.type));
end

if isfield(input_scenario, "initial_position_I")
    S.initial_position_I = input_scenario.initial_position_I(:);
end
if isfield(input_scenario, "initial_velocity_I")
    S.initial_velocity_I = input_scenario.initial_velocity_I(:);
end
if isfield(input_scenario, "initial_quaternion_BI")
    S.initial_quaternion_BI = input_scenario.initial_quaternion_BI(:);
end
if isfield(input_scenario, "initial_omega_B")
    S.initial_omega_B = input_scenario.initial_omega_B(:);
end
if isfield(input_scenario, "initial_actuators")
    S.initial_actuators = input_scenario.initial_actuators(:);
end

if P.metrics.evaluate_after_entry
    S.metric_start_time = S.entry_time;
else
    S.metric_start_time = 0;
end

require_positive(P.simulation.controller_sample_time, ...
    "P.simulation.controller_sample_time");
require_positive(P.simulation.integration_substeps, ...
    "P.simulation.integration_substeps");
if P.simulation.integration_substeps ~= ...
        round(P.simulation.integration_substeps)
    error('simulate_COMET:InvalidIntegrationSubsteps', ...
        "P.simulation.integration_substeps must be an integer.");
end
require_positive(S.final_time, "scenario final time");

end


function entry_time = choose_entry_time(P, input_scenario, S)

if ~S.start_from_rest || ~P.trajectory.use_acceleration_limited_entry
    entry_time = 0;
    return
end

if isfield(input_scenario, "entry_time")
    entry_time = input_scenario.entry_time;
    if entry_time < 0
        error('simulate_COMET:InvalidEntryTime', ...
            "scenario.entry_time cannot be negative.");
    end
    return
end

default_entry_acceleration = max(1.5, min( ...
    P.controller.velocity.acceleration_limit(1), max(S.amax, 1.5)));
entry_acceleration = get_option(input_scenario, ...
    "entry_acceleration", default_entry_acceleration);
require_positive(entry_acceleration, "scenario.entry_acceleration");

% For the quintic smoothstep v/vmax = 10 xi^3 - 15 xi^4 + 6 xi^5, the
% maximum tangential acceleration is 1.875*vmax/entry_time.
entry_time = 1.875*S.vmax/entry_acceleration;

% Avoid an entry ramp that would consume the whole requested trajectory.
maximum_entry_time = S.target_phase / max(S.phase_rate, eps);
entry_time = min(entry_time, maximum_entry_time);

end


function final_time = choose_final_time(input_scenario, S)

if isfield(input_scenario, "duration")
    final_time = input_scenario.duration;
    return
end

% The smooth entry covers 0.5*phase_rate*entry_time radians. Therefore
% target_phase is reached at target/phase_rate + entry_time/2.
final_time = S.target_phase/S.phase_rate + 0.5*S.entry_time;

end


function value = get_option(S, field_name, default_value)

if isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end

end


function require_positive(value, label)

if ~isscalar(value) || ~isfinite(value) || value <= 0
    error('simulate_COMET:InvalidPositiveValue', ...
        "%s must be a finite positive scalar.", label);
end

end


%% Initialization and logging

function x = initial_state(P, S, I)

x = zeros(I.n_states, 1);
ref0 = trajectory_reference(P, S, 0);

x(I.position) = ref0.position_I;
x(I.velocity) = ref0.velocity_I;
x(I.quaternion) = [1; 0; 0; 0];
x(I.omega) = zeros(3,1);
x(I.actuators) = [ ...
    P.physical.hover_thrust_per_rotor;
    P.physical.hover_thrust_per_rotor;
    0; 0; 0; 0];

if isfield(S, "initial_position_I")
    assert(numel(S.initial_position_I) == 3, ...
        "initial_position_I must contain three elements.");
    x(I.position) = S.initial_position_I;
end
if isfield(S, "initial_velocity_I")
    assert(numel(S.initial_velocity_I) == 3, ...
        "initial_velocity_I must contain three elements.");
    x(I.velocity) = S.initial_velocity_I;
end
if isfield(S, "initial_quaternion_BI")
    assert(numel(S.initial_quaternion_BI) == 4, ...
        "initial_quaternion_BI must contain four elements.");
    x(I.quaternion) = normalize_quaternion(S.initial_quaternion_BI);
end
if isfield(S, "initial_omega_B")
    assert(numel(S.initial_omega_B) == 3, ...
        "initial_omega_B must contain three elements.");
    x(I.omega) = S.initial_omega_B;
end
if isfield(S, "initial_actuators")
    assert(numel(S.initial_actuators) == I.n_actuators, ...
        "initial_actuators must contain six elements.");
    x(I.actuators) = S.initial_actuators;
end

end


function memory = initial_controller_memory()

memory.velocity_integral = zeros(3,1);
memory.rate_integral = zeros(3,1);
memory.velocity_derivative_filtered = zeros(3,1);
memory.rate_derivative_filtered = zeros(3,1);
memory.previous_velocity_error = zeros(3,1);
memory.previous_rate_error = zeros(3,1);
memory.previous_desired_rotation = eye(3);
memory.first_sample = true;

end


function log = preallocate_log(N, n_actuators)

log.state = nan(19, N);

log.reference.position_I = nan(3, N);
log.reference.velocity_I = nan(3, N);
log.reference.acceleration_I = nan(3, N);
log.reference.yaw = nan(1, N);
log.reference.phase = nan(1, N);

log.error.position_I = nan(3, N);
log.error.velocity_I = nan(3, N);
log.error.attitude = nan(3, N);
log.error.angular_rate_B = nan(3, N);

log.control.velocity_command_I = nan(3, N);
log.control.acceleration_command_I = nan(3, N);
log.control.velocity_error_derivative_raw = nan(3, N);
log.control.velocity_error_derivative_filtered = nan(3, N);
log.control.force_desired_I = nan(3, N);
log.control.quaternion_desired_BI = nan(4, N);
log.control.omega_desired_B = nan(3, N);
log.control.moment_desired_B = nan(3, N);
log.control.rate_error_derivative_raw = nan(3, N);
log.control.rate_error_derivative_filtered = nan(3, N);
log.control.wrench_desired_B = nan(6, N);

log.actuator.command_raw = nan(n_actuators, N);
log.actuator.command_limited = nan(n_actuators, N);
log.actuator.actual = nan(n_actuators, N);
log.actuator.rate_raw = nan(n_actuators, N);
log.actuator.rate_limited = nan(n_actuators, N);

log.saturation.position = false(n_actuators, N);
log.saturation.rate = false(n_actuators, N);

log.allocation.cyclic_effectiveness = nan(4, N);
log.allocation.reciprocal_condition = nan(1, N);
log.allocation.used_pseudoinverse = false(1, N);
log.allocation.moment_residual = nan(2, N);

log.wrench.rotor_B = nan(6, N);
log.wrench.rotor_force_I = nan(3, N);
log.wrench.environment_force_I = nan(3, N);
log.wrench.environment_moment_B = nan(3, N);
log.wrench.residual_B = nan(6, N);

log.rotor.force_up_B = nan(3, N);
log.rotor.force_dw_B = nan(3, N);
log.rotor.direction_up_B = nan(3, N);
log.rotor.direction_dw_B = nan(3, N);

log.environment.drag_force_I = nan(3, N);
log.environment.gust_force_I = nan(3, N);
log.environment.drag_moment_B = nan(3, N);
log.environment.disturbance_moment_B = nan(3, N);

log.power.estimated_total_W = nan(1, N);

end


function log = store_sample(log, k, x, ref, control, allocation, ...
    rate_status, rotor_wrench, rotor_detail, environment, ...
    environment_detail, P, I)

q = normalize_quaternion(x(I.quaternion));
R = quaternion_to_rotation(q);

log.state(:,k) = x;
log.state(I.quaternion,k) = q;

log.reference.position_I(:,k) = ref.position_I;
log.reference.velocity_I(:,k) = ref.velocity_I;
log.reference.acceleration_I(:,k) = ref.acceleration_I;
log.reference.yaw(k) = ref.yaw;
log.reference.phase(k) = ref.phase;

% Tracking metrics must use the current physical state, not the delayed
% state internally seen by the controller. Otherwise increasing the delay
% would contaminate the reported error with a measurement-age artefact.
log.error.position_I(:,k) = ref.position_I - x(I.position);
log.error.velocity_I(:,k) = ref.velocity_I - x(I.velocity);
log.error.attitude(:,k) = control.attitude_error;
log.error.angular_rate_B(:,k) = control.rate_error_B;

log.control.velocity_command_I(:,k) = control.velocity_command_I;
log.control.acceleration_command_I(:,k) = control.acceleration_command_I;
log.control.velocity_error_derivative_raw(:,k) = ...
    control.velocity_error_derivative_raw;
log.control.velocity_error_derivative_filtered(:,k) = ...
    control.velocity_error_derivative_filtered;
log.control.force_desired_I(:,k) = control.force_desired_I;
log.control.quaternion_desired_BI(:,k) = control.quaternion_desired_BI;
log.control.omega_desired_B(:,k) = control.omega_desired_B;
log.control.moment_desired_B(:,k) = control.moment_desired_B;
log.control.rate_error_derivative_raw(:,k) = ...
    control.rate_error_derivative_raw;
log.control.rate_error_derivative_filtered(:,k) = ...
    control.rate_error_derivative_filtered;
log.control.wrench_desired_B(:,k) = control.wrench_desired_B;

log.actuator.command_raw(:,k) = allocation.command_raw;
log.actuator.command_limited(:,k) = allocation.command_limited;
log.actuator.actual(:,k) = bounded_actuator_state(P, x(I.actuators));
log.actuator.rate_raw(:,k) = rate_status.raw;
log.actuator.rate_limited(:,k) = rate_status.limited;
log.saturation.position(:,k) = allocation.position_saturated;
log.saturation.rate(:,k) = rate_status.saturated;
log.allocation.cyclic_effectiveness(:,k) = ...
    allocation.cyclic_effectiveness_matrix(:);
log.allocation.reciprocal_condition(k) = ...
    allocation.reciprocal_condition;
log.allocation.used_pseudoinverse(k) = ...
    allocation.used_pseudoinverse;
log.allocation.moment_residual(:,k) = ...
    allocation.moment_residual;

log.wrench.rotor_B(:,k) = rotor_wrench;
log.wrench.rotor_force_I(:,k) = R * rotor_wrench(1:3);
log.wrench.environment_force_I(:,k) = environment.force_I;
log.wrench.environment_moment_B(:,k) = environment.moment_B;
log.wrench.residual_B(:,k) = rotor_wrench - control.wrench_desired_B;

log.rotor.force_up_B(:,k) = rotor_detail.force_up_B;
log.rotor.force_dw_B(:,k) = rotor_detail.force_dw_B;
log.rotor.direction_up_B(:,k) = rotor_detail.direction_up_B;
log.rotor.direction_dw_B(:,k) = rotor_detail.direction_dw_B;

log.environment.drag_force_I(:,k) = environment_detail.drag_force_I;
log.environment.gust_force_I(:,k) = environment_detail.gust_force_I;
log.environment.drag_moment_B(:,k) = environment_detail.drag_moment_B;
log.environment.disturbance_moment_B(:,k) = ...
    environment_detail.disturbance_moment_B;

log.power.estimated_total_W(k) = estimate_power(P, x(I.actuators));

end


%% Reference trajectory

function ref = trajectory_reference(~, S, t)

[theta, theta_dot, theta_ddot] = phase_profile(S, t);

switch S.type
    case "hover"
        position = [0; 0; S.height];
        velocity = zeros(3,1);
        acceleration = zeros(3,1);

    case "circle"
        c = cos(theta);
        s = sin(theta);

        position = [S.radius*c; S.radius*s; S.height];
        velocity = S.radius * theta_dot * [-s; c; 0];
        acceleration = S.radius * ( ...
            theta_ddot * [-s; c; 0] - ...
            theta_dot^2 * [c; s; 0]);

    case "eight"
        position = [ ...
            S.A*sin(theta);
            S.B*sin(2*theta);
            S.height];

        velocity = [ ...
            S.A*cos(theta)*theta_dot;
            2*S.B*cos(2*theta)*theta_dot;
            0];

        acceleration = [ ...
            S.A*(-sin(theta)*theta_dot^2 + cos(theta)*theta_ddot);
            2*S.B*(-2*sin(2*theta)*theta_dot^2 + ...
                    cos(2*theta)*theta_ddot);
            0];

    otherwise
        error('simulate_COMET:InternalTrajectoryError', ...
            "Unsupported trajectory type.");
end

ref.position_I = position;
ref.velocity_I = velocity;
ref.acceleration_I = acceleration;
ref.yaw = S.yaw;
ref.phase = theta;
ref.phase_rate = theta_dot;
ref.phase_acceleration = theta_ddot;

end


function [theta, theta_dot, theta_ddot] = phase_profile(S, t)

if S.type == "hover"
    theta = 0;
    theta_dot = 0;
    theta_ddot = 0;
    return
end

omega = S.phase_rate;
T = S.entry_time;

if T > 0 && t < T
    xi = t/T;
    % Quintic velocity scaling. Velocity, acceleration and jerk match the
    % constant-speed segment continuously at the end of the entry phase.
    velocity_scale = 10*xi^3 - 15*xi^4 + 6*xi^5;
    theta = omega*T*(2.5*xi^4 - 3*xi^5 + xi^6);
    theta_dot = omega*velocity_scale;
    theta_ddot = (omega/T)*(30*xi^2 - 60*xi^3 + 30*xi^4);
else
    theta_at_entry = 0.5*omega*T;
    theta = theta_at_entry + omega*(t - T);
    theta_dot = omega;
    theta_ddot = 0;
end

end


%% Cascaded controller

function [control, sample] = controller_sample(P, x, ref, memory, dt, I)

position = x(I.position);
velocity = x(I.velocity);
q = normalize_quaternion(x(I.quaternion));
omega = x(I.omega);
R = quaternion_to_rotation(q);

position_error = ref.position_I - position;
velocity_command_raw = ref.velocity_I + ...
    P.controller.position.Kp * position_error;
velocity_command = component_limit(velocity_command_raw, ...
    P.controller.position.velocity_command_limit);

velocity_error = velocity_command - velocity;
if memory.first_sample
    velocity_error_derivative_raw = zeros(3,1);
    velocity_error_derivative_filtered = zeros(3,1);
else
    velocity_error_derivative_raw = ...
        (velocity_error - memory.previous_velocity_error)/dt;
    alpha_velocity = exp(-dt / ...
        P.controller.velocity.derivative_filter_time_constant);
    velocity_error_derivative_filtered = ...
        alpha_velocity*memory.velocity_derivative_filtered + ...
        (1-alpha_velocity)*velocity_error_derivative_raw;
end

acceleration_command_raw = ref.acceleration_I + ...
    P.controller.velocity.Kp * velocity_error + ...
    P.controller.velocity.Ki * memory.velocity_integral + ...
    P.controller.velocity.Kd * velocity_error_derivative_filtered;
acceleration_command = component_limit(acceleration_command_raw, ...
    P.controller.velocity.acceleration_limit);

gravity_I = P.physical.g * P.frames.gravity_axis_I;
force_desired_I = P.physical.mass * ...
    (acceleration_command - gravity_I);

[R_desired, q_desired] = desired_attitude_from_force( ...
    force_desired_I, ref.yaw);

q_error = quaternion_multiply(quaternion_conjugate(q), q_desired);
if q_error(1) < 0
    q_error = -q_error;
end
attitude_error = 2*q_error(2:4);

if memory.first_sample
    omega_feedforward_current_B = zeros(3,1);
else
    rotation_increment = memory.previous_desired_rotation.' * R_desired;
    omega_desired_desired_B = rotation_log_vector(rotation_increment)/dt;
    omega_feedforward_current_B = ...
        R.' * R_desired * omega_desired_desired_B;
end

omega_desired_raw = omega_feedforward_current_B + ...
    P.controller.attitude.Kp * attitude_error;
omega_desired = component_limit(omega_desired_raw, ...
    P.controller.attitude.omega_command_limit);

rate_error = omega_desired - omega;
if memory.first_sample
    rate_error_derivative_raw = zeros(3,1);
    rate_error_derivative_filtered = zeros(3,1);
else
    rate_error_derivative_raw = ...
        (rate_error - memory.previous_rate_error)/dt;
    alpha_rate = exp(-dt / ...
        P.controller.rate.derivative_filter_time_constant);
    rate_error_derivative_filtered = ...
        alpha_rate*memory.rate_derivative_filtered + ...
        (1-alpha_rate)*rate_error_derivative_raw;
end

% cross(omega,J*omega) compensates the rigid-body gyroscopic term. The
% remaining PID contribution creates the requested angular acceleration.
moment_desired = cross(omega, P.physical.J*omega) + ...
    P.controller.rate.Kp * rate_error + ...
    P.controller.rate.Ki * memory.rate_integral + ...
    P.controller.rate.Kd * rate_error_derivative_filtered;

force_desired_B = R.' * force_desired_I;

control.position_error_I = position_error;
control.velocity_error_I = velocity_error;
control.velocity_command_I = velocity_command;
control.acceleration_command_I = acceleration_command;
control.velocity_error_derivative_raw = velocity_error_derivative_raw;
control.velocity_error_derivative_filtered = ...
    velocity_error_derivative_filtered;
control.force_desired_I = force_desired_I;
control.total_thrust_desired = norm(force_desired_I);
control.rotation_desired_BI = R_desired;
control.quaternion_desired_BI = q_desired;
control.attitude_error = attitude_error;
control.omega_desired_B = omega_desired;
control.rate_error_B = rate_error;
control.moment_desired_B = moment_desired;
control.rate_error_derivative_raw = rate_error_derivative_raw;
control.rate_error_derivative_filtered = rate_error_derivative_filtered;
control.wrench_desired_B = [force_desired_B; moment_desired];

sample = memory;
sample.current_velocity_error = velocity_error;
sample.current_rate_error = rate_error;
sample.current_desired_rotation = R_desired;
sample.current_velocity_derivative_filtered = ...
    velocity_error_derivative_filtered;
sample.current_rate_derivative_filtered = ...
    rate_error_derivative_filtered;

end


function [R_desired, q_desired] = desired_attitude_from_force(force_I, yaw)

force_norm = norm(force_I);
if force_norm < 1e-9
    b3_desired = [0; 0; 1];
else
    b3_desired = force_I/force_norm;
end

b1_heading = [cos(yaw); sin(yaw); 0];
b2_desired = cross(b3_desired, b1_heading);

if norm(b2_desired) < 1e-8
    alternate_heading = [-sin(yaw); cos(yaw); 0];
    b2_desired = cross(b3_desired, alternate_heading);
end

b2_desired = b2_desired/norm(b2_desired);
b1_desired = cross(b2_desired, b3_desired);
b1_desired = b1_desired/norm(b1_desired);

R_desired = [b1_desired, b2_desired, b3_desired];
q_desired = rotation_to_quaternion(R_desired);

end


function memory = update_controller_memory(P, sample, control, ...
    saturation_active, dt)

if saturation_active
    integration_scale = max(0, 1 - P.controller.anti_windup_gain);
else
    integration_scale = 1;
end

velocity_integral = sample.velocity_integral + ...
    integration_scale*dt*control.velocity_error_I;
rate_integral = sample.rate_integral + ...
    integration_scale*dt*control.rate_error_B;

memory.velocity_integral = component_limit(velocity_integral, ...
    P.controller.velocity.integral_limit);
memory.rate_integral = component_limit(rate_integral, ...
    P.controller.rate.integral_limit);
memory.velocity_derivative_filtered = ...
    sample.current_velocity_derivative_filtered;
memory.rate_derivative_filtered = ...
    sample.current_rate_derivative_filtered;
memory.previous_velocity_error = sample.current_velocity_error;
memory.previous_rate_error = sample.current_rate_error;
memory.previous_desired_rotation = sample.current_desired_rotation;
memory.first_sample = false;

end


%% Control allocation and actuator dynamics

function allocation = allocate_actuators(P, C, control)

T_total = control.total_thrust_desired;
moment = control.moment_desired_B;

k_yaw = P.mixer.yaw_moment_per_thrust_difference;
if abs(k_yaw) < 1e-9
    thrust_difference = 0;
else
    thrust_difference = moment(3)/k_yaw;
end

% Positive yaw moment is defined as a larger lower-rotor reaction torque.
T_up_raw = 0.5*(T_total - thrust_difference);
T_dw_raw = 0.5*(T_total + thrust_difference);

% Evaluate cyclic authority at the thrusts that the physical actuators can
% actually produce. The high-level controller remains identical for both
% configurations; only this actuator-to-wrench map changes.
T_up_for_map = min(max(T_up_raw, P.actuators.thrust_min), ...
    P.actuators.thrust_max_per_rotor);
T_dw_for_map = min(max(T_dw_raw, P.actuators.thrust_min), ...
    P.actuators.thrust_max_per_rotor);

[cyclic_map, baseline_moment] = cyclic_effectiveness_matrix( ...
    P, C, T_up_for_map, T_dw_for_map);
required_cyclic_moment = moment(1:2) - baseline_moment;
reciprocal_condition = rcond(cyclic_map);

use_pseudoinverse = reciprocal_condition < ...
    P.allocation.minimum_reciprocal_condition;
if use_pseudoinverse
    if ~P.allocation.use_pseudoinverse_if_ill_conditioned
        error('simulate_COMET:SingularAllocation', ...
            "The cyclic-effectiveness matrix is singular or ill-conditioned.");
    end
    common_cyclic_raw = pinv(cyclic_map) * required_cyclic_moment;
else
    common_cyclic_raw = cyclic_map \ required_cyclic_moment;
end

aileron_common_raw = common_cyclic_raw(1);
elevator_common_raw = common_cyclic_raw(2);

cyclic_raw = [ ...
    aileron_common_raw;
    elevator_common_raw;
    aileron_common_raw;
    elevator_common_raw] .* C.cyclic_mask;

command_raw = [T_up_raw; T_dw_raw; cyclic_raw];
lower = [ ...
    P.actuators.thrust_min;
    P.actuators.thrust_min;
    repmat(P.actuators.cyclic_min, 4, 1)];
upper = [ ...
    P.actuators.thrust_max_per_rotor;
    P.actuators.thrust_max_per_rotor;
    repmat(P.actuators.cyclic_max, 4, 1)];

command_limited = min(max(command_raw, lower), upper);
command_limited(3:6) = command_limited(3:6).*C.cyclic_mask;

active_channels = [true; true; logical(C.cyclic_mask)];
position_saturated = active_channels & ...
    abs(command_raw - command_limited) > ...
    P.actuators.saturation_tolerance;

allocation.command_raw = command_raw;
allocation.command_limited = command_limited;
allocation.position_saturated = position_saturated;
allocation.common_cyclic_raw = [aileron_common_raw; elevator_common_raw];
allocation.cyclic_effectiveness_matrix = cyclic_map;
allocation.reciprocal_condition = reciprocal_condition;
allocation.used_pseudoinverse = use_pseudoinverse;
allocation.moment_residual = ...
    cyclic_map*common_cyclic_raw - required_cyclic_moment;

end


function [A, baseline_moment] = cyclic_effectiveness_matrix( ...
    P, C, T_up, T_dw)
% Local map from common aileron/elevator commands to roll/pitch moments.
%
% The derivative is computed from the complete nonlinear rotor-wrench model,
% so it automatically includes configuration masks, rotor interaction,
% flapping/hub coefficients, thrust-vector tilt and hub moment arms.

h = P.allocation.cyclic_linearization_step;
u0 = [T_up; T_dw; 0; 0; 0; 0];
[wrench_0, ~] = rotor_wrench_body(P, C, u0);
baseline_moment = wrench_0(4:5);

A = zeros(2,2);
for j = 1:2
    cyclic_direction = zeros(4,1);
    cyclic_direction(j) = 1;
    cyclic_direction(j+2) = 1;

    u_plus = [T_up; T_dw;  h*cyclic_direction];
    u_minus = [T_up; T_dw; -h*cyclic_direction];
    [wrench_plus, ~] = rotor_wrench_body(P, C, u_plus);
    [wrench_minus, ~] = rotor_wrench_body(P, C, u_minus);

    A(:,j) = (wrench_plus(4:5) - wrench_minus(4:5))/(2*h);
end

end


function status = actuator_rate(P, C, actuator_state, command)

actual = bounded_actuator_state(P, actuator_state);

raw = zeros(6,1);
raw(1:2) = (command(1:2) - actual(1:2)) / ...
    P.actuators.motor_time_constant;
raw(3:6) = (command(3:6) - actual(3:6)) / ...
    P.actuators.servo_time_constant;

limits = [ ...
    P.actuators.thrust_rate_max;
    P.actuators.thrust_rate_max;
    repmat(P.actuators.cyclic_rate_max, 4, 1)];

if P.actuators.apply_rate_limits
    limited = min(max(raw, -limits), limits);
else
    limited = raw;
end

% Do not integrate farther outside a hard position bound.
at_lower_thrust = actual(1:2) <= P.actuators.thrust_min;
at_upper_thrust = actual(1:2) >= P.actuators.thrust_max_per_rotor;
limited(1:2) = block_outward_rate(limited(1:2), ...
    at_lower_thrust, at_upper_thrust);

at_lower_cyclic = actual(3:6) <= P.actuators.cyclic_min;
at_upper_cyclic = actual(3:6) >= P.actuators.cyclic_max;
limited(3:6) = block_outward_rate(limited(3:6), ...
    at_lower_cyclic, at_upper_cyclic);
limited(3:6) = limited(3:6).*C.cyclic_mask;

active_channels = [true; true; logical(C.cyclic_mask)];
status.raw = raw;
status.limited = limited;
status.saturated = active_channels & ...
    abs(raw - limited) > P.actuators.saturation_tolerance;

end


function rate = block_outward_rate(rate, at_lower, at_upper)

rate(at_lower & rate < 0) = 0;
rate(at_upper & rate > 0) = 0;

end


function actual = bounded_actuator_state(P, actuator_state)

actual = actuator_state(:);
actual(1:2) = min(max(actual(1:2), P.actuators.thrust_min), ...
    P.actuators.thrust_max_per_rotor);
actual(3:6) = min(max(actual(3:6), P.actuators.cyclic_min), ...
    P.actuators.cyclic_max);

end


%% Rotor and environmental wrench

function [wrench_B, detail] = rotor_wrench_body(P, C, actuator_state)

u = bounded_actuator_state(P, actuator_state);
T_up = u(1);
T_dw = u(2);
ail_up = u(3)*C.cyclic_mask(1);
ele_up = u(4)*C.cyclic_mask(2);
ail_dw = u(5)*C.cyclic_mask(3);
ele_dw = u(6)*C.cyclic_mask(4);

% The normalized servo positions are converted into effective aerodynamic
% cyclic deflections through one common calibration factor. Applying the
% factor here keeps the nonlinear plant and the configuration-dependent
% allocation map consistent, because both use rotor_wrench_body().
K_cyclic = P.aerodynamics.cyclic_effectiveness_scale;
ail_up_effective = K_cyclic*ail_up;
ele_up_effective = K_cyclic*ele_up;
ail_dw_effective = K_cyclic*ail_dw;
ele_dw_effective = K_cyclic*ele_dw;

direction_up = rotor_thrust_direction( ...
    P.aerodynamics.A_up*ele_up_effective, ...
    P.aerodynamics.B_up*ail_up_effective);
direction_dw = rotor_thrust_direction( ...
    P.aerodynamics.A_dw*ele_dw_effective, ...
    P.aerodynamics.B_dw*ail_dw_effective);

force_up = T_up*direction_up;
force_dw = T_dw*direction_dw;
force_total = force_up + force_dw;

moment_arm_up = cross(P.geometry.hub_position_up_B, force_up);
moment_arm_dw = cross(P.geometry.hub_position_dw_B, force_dw);

if C.is_dual
    interaction = P.aerodynamics.dual_interaction_factor;
else
    interaction = 1.0;
end

% Residual flapping/hub moments. The geometric r x F terms above are not
% included in these coefficients and are therefore added only once.
moment_flap_up = interaction*T_up*[ ...
    P.aerodynamics.cyclic_arm_roll_up*ail_up_effective;
    P.aerodynamics.cyclic_arm_pitch_up*ele_up_effective;
    0];
moment_flap_dw = interaction*T_dw*[ ...
    P.aerodynamics.cyclic_arm_roll_dw*ail_dw_effective;
    P.aerodynamics.cyclic_arm_pitch_dw*ele_dw_effective;
    0];

moment_reaction = [0; 0; ...
    P.aerodynamics.yaw_torque_per_thrust*(T_dw - T_up)];

moment_total = moment_arm_up + moment_arm_dw + ...
    moment_flap_up + moment_flap_dw + moment_reaction;

wrench_B = [force_total; moment_total];

detail.force_up_B = force_up;
detail.force_dw_B = force_dw;
detail.direction_up_B = direction_up;
detail.direction_dw_B = direction_dw;
detail.moment_arm_up_B = moment_arm_up;
detail.moment_arm_dw_B = moment_arm_dw;
detail.moment_flap_up_B = moment_flap_up;
detail.moment_flap_dw_B = moment_flap_dw;
detail.moment_reaction_B = moment_reaction;

end


function direction = rotor_thrust_direction(elevator_angle, aileron_angle)

% Equivalent to R_x(aileron_angle)*R_y(elevator_angle)*e3. Positive
% elevator tilts thrust toward +x_B; positive aileron tilts it toward -y_B.
direction = [ ...
    sin(elevator_angle);
    -sin(aileron_angle)*cos(elevator_angle);
    cos(aileron_angle)*cos(elevator_angle)];
direction = direction/norm(direction);

end


function [environment, detail] = environmental_wrench(P, t, x, I)

q = normalize_quaternion(x(I.quaternion));
R = quaternion_to_rotation(q);
velocity_B = R.'*x(I.velocity);
omega_B = x(I.omega);

drag_force_B = -P.aerodynamics.linear_drag_B*velocity_B - ...
    P.aerodynamics.quadratic_drag_B*(abs(velocity_B).*velocity_B);
drag_force_I = R*drag_force_B;
drag_moment_B = -P.aerodynamics.angular_damping_B*omega_B;

gust_force_I = zeros(3,1);
disturbance_moment_B = zeros(3,1);

if P.disturbance.enabled
    disturbance_type = lower(string(P.disturbance.type));

    if disturbance_type == "gust" || disturbance_type == "both"
        gust_force_I = P.disturbance.gust.force_I * ...
            one_cosine_pulse(t, P.disturbance.gust.start_time, ...
            P.disturbance.gust.duration);
    end

    if disturbance_type == "torque" || disturbance_type == "both"
        active = t >= P.disturbance.torque.start_time && ...
            t <= P.disturbance.torque.start_time + ...
            P.disturbance.torque.duration;
        if active
            disturbance_moment_B = P.disturbance.torque.value_B;
        end
    end
end

environment.force_I = drag_force_I + gust_force_I;
environment.moment_B = drag_moment_B + disturbance_moment_B;

detail.drag_force_I = drag_force_I;
detail.gust_force_I = gust_force_I;
detail.drag_moment_B = drag_moment_B;
detail.disturbance_moment_B = disturbance_moment_B;

end


function amplitude = one_cosine_pulse(t, start_time, duration)

if duration <= 0 || t < start_time || t > start_time + duration
    amplitude = 0;
    return
end

xi = (t - start_time)/duration;
amplitude = 0.5*(1 - cos(2*pi*xi));

end


function power_W = estimate_power(P, actuator_state)

u = bounded_actuator_state(P, actuator_state);
T_reference = 0.5*P.physical.mass_endurance_standard*P.physical.g;

% Reduced induced-power model P proportional to T^(3/2), normalized so that
% equal reference thrust reproduces the reported 163.2 W hover measurement.
power_W = 0.5*P.propulsion.hover_power_standard * ( ...
    (u(1)/T_reference)^(3/2) + ...
    (u(2)/T_reference)^(3/2));

end


%% Continuous dynamics and RK4

function x_dot = state_derivative(t, x, actuator_command, P, C, I)

q = normalize_quaternion(x(I.quaternion));
R = quaternion_to_rotation(q);
omega = x(I.omega);

[rotor_wrench, ~] = rotor_wrench_body(P, C, x(I.actuators));
[environment, ~] = environmental_wrench(P, t, x, I);

gravity_I = P.physical.g*P.frames.gravity_axis_I;
force_rotor_I = R*rotor_wrench(1:3);
acceleration_I = gravity_I + ...
    (force_rotor_I + environment.force_I)/P.physical.mass;

moment_total_B = rotor_wrench(4:6) + environment.moment_B;
omega_dot = P.physical.J \ ( ...
    moment_total_B - cross(omega, P.physical.J*omega));

q_dot = 0.5*quaternion_multiply(q, [0; omega]);
actuator_status = actuator_rate(P, C, x(I.actuators), ...
    actuator_command);

x_dot = zeros(I.n_states,1);
x_dot(I.position) = x(I.velocity);
x_dot(I.velocity) = acceleration_I;
x_dot(I.quaternion) = q_dot;
x_dot(I.omega) = omega_dot;
x_dot(I.actuators) = actuator_status.limited;

end


function x_next = rk4_step(dynamics, t, x, h)

k1 = dynamics(t,         x);
k2 = dynamics(t + h/2,   x + (h/2)*k1);
k3 = dynamics(t + h/2,   x + (h/2)*k2);
k4 = dynamics(t + h,     x + h*k3);

x_next = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

end


function [x, ok] = project_state(P, C, x, I)

ok = all(isfinite(x));
if ~ok
    return
end

q_norm = norm(x(I.quaternion));
if q_norm < 1e-12
    ok = false;
    return
end

if P.simulation.normalize_quaternion_each_step
    x(I.quaternion) = x(I.quaternion)/q_norm;
end

x(I.actuators) = bounded_actuator_state(P, x(I.actuators));
x(I.actuators(3:6)) = ...
    x(I.actuators(3:6)).*C.cyclic_mask;

end


function [ok, reason] = validate_dynamic_state(P, x, I)

ok = false;
reason = "unknown_state_failure";

if any(~isfinite(x))
    reason = "non_finite_state";
    return
end
if norm(x(I.quaternion)) < 1e-12
    reason = "invalid_quaternion";
    return
end
if norm(x(I.position)) > P.simulation.max_position_norm
    reason = "position_divergence";
    return
end
if norm(x(I.velocity)) > P.simulation.max_velocity_norm
    reason = "velocity_divergence";
    return
end
if norm(x(I.omega)) > P.simulation.max_angular_rate
    reason = "angular_rate_divergence";
    return
end

ok = true;
reason = "state_valid";

end


%% Output and metrics

function out = assemble_output(P, C, S, t, log, last_index, ...
    termination_reason, numerical_divergence, stopped_for_lost_control)

fields = fieldnames(log);
for i = 1:numel(fields)
    log.(fields{i}) = trim_nested(log.(fields{i}), last_index);
end
t = t(1:last_index);

out.time = t;
out.configuration = C.name;
out.cyclic_mask = C.cyclic_mask;
out.scenario = S;
out.parameter_version = P.meta.version;
out.parameters = P;

out.state.position_I = log.state(1:3,:);
out.state.velocity_I = log.state(4:6,:);
out.state.quaternion_BI = log.state(7:10,:);
out.state.omega_B = log.state(11:13,:);
out.state.actuators = log.state(14:19,:);

out.reference = log.reference;
out.error = log.error;
out.control = log.control;
out.actuator = log.actuator;
out.saturation = log.saturation;
out.allocation = log.allocation;
out.wrench = log.wrench;
out.rotor = log.rotor;
out.environment = log.environment;
out.power = log.power;

out.status.termination_reason = termination_reason;
out.status.numerical_divergence = numerical_divergence;
out.status.stopped_for_lost_control = stopped_for_lost_control;
out.status.simulated_final_time = t(end);
out.status.requested_final_time = S.final_time;

out.metrics = compute_metrics(P, C, S, t, log, numerical_divergence);
out.status.lost_control = out.metrics.lost_control;
out.status.completed = ~numerical_divergence && ...
    ~out.metrics.lost_control;

out.notes.integrator = ...
    "Fixed-step classical RK4 with zero-order-held controller output.";
out.notes.controller_sample_time_s = ...
    P.simulation.controller_sample_time;
out.notes.integration_substeps = P.simulation.integration_substeps;
out.notes.integration_step_s = ...
    P.simulation.controller_sample_time/P.simulation.integration_substeps;
out.notes.feedback_delay_requested_s = ...
    P.feedback.state_delay_requested;
out.notes.feedback_delay_applied_s = P.feedback.state_delay_applied;
out.notes.feedback_delay_samples = P.feedback.state_delay_samples;
out.notes.feedback_delay_floor_samples = ...
    P.feedback.state_delay_floor_samples;
out.notes.feedback_delay_fraction = P.feedback.state_delay_fraction;
out.notes.derivative_filter = ...
    "First-order low-pass filters are applied to discrete velocity- and rate-error derivatives.";
out.notes.power_model = ...
    "Estimated from P proportional to T^(3/2); not a measured electrical model.";
out.notes.cyclic_model = ...
    "cyclic_arm fields are residual flapping coefficients; r_hub cross F is added separately.";
out.notes.control_allocation = ...
    "A configuration-dependent local cyclic-effectiveness matrix is computed from the nonlinear rotor-wrench model at every controller sample.";

end


function value = trim_nested(value, last_index)

if isstruct(value)
    names = fieldnames(value);
    for j = 1:numel(names)
        value.(names{j}) = trim_nested(value.(names{j}), last_index);
    end
elseif isnumeric(value) || islogical(value)
    if ~isscalar(value) && size(value,2) >= last_index
        value = value(:,1:last_index);
    end
end

end


function metrics = compute_metrics(P, C, S, t, log, numerical_divergence)

evaluation = t >= S.metric_start_time;
if ~any(evaluation)
    evaluation(:) = true;
end

position_error = log.error.position_I(:,evaluation);
velocity_error = log.error.velocity_I(:,evaluation);
position_error_norm = vecnorm(position_error, 2, 1);
velocity_error_norm = vecnorm(velocity_error, 2, 1);

metrics.metric_start_time = S.metric_start_time;
metrics.samples_evaluated = nnz(evaluation);
metrics.rmse_position = sqrt(mean(sum(position_error.^2,1)));
metrics.rmse_position_axis = sqrt(mean(position_error.^2,2));
metrics.max_position_error = max(position_error_norm);
metrics.rmse_velocity = sqrt(mean(sum(velocity_error.^2,1)));
metrics.max_velocity_error = max(velocity_error_norm);

position_sat = log.saturation.position(:,evaluation);
rate_sat = log.saturation.rate(:,evaluation);

metrics.thrust_position_saturation_fraction = ...
    mean(any(position_sat(1:2,:),1));
metrics.cyclic_position_saturation_fraction = ...
    mean(any(position_sat(3:6,:),1));
metrics.actuator_rate_saturation_fraction = ...
    mean(any(rate_sat,1));
metrics.any_saturation_fraction = ...
    mean(any(position_sat | rate_sat,1));

metrics.minimum_allocation_reciprocal_condition = ...
    min(log.allocation.reciprocal_condition(evaluation));
metrics.allocation_pseudoinverse_fraction = ...
    mean(log.allocation.used_pseudoinverse(evaluation));
metrics.allocation_moment_residual_rms = sqrt(mean( ...
    sum(log.allocation.moment_residual(:,evaluation).^2,1)));

roll_rows = [3,5];
roll_rows = roll_rows(logical(C.cyclic_mask([1,3])));
pitch_rows = [4,6];
pitch_rows = pitch_rows(logical(C.cyclic_mask([2,4])));

metrics.roll_command_range = [ ...
    min(log.actuator.command_limited(roll_rows,evaluation), [], "all"), ...
    max(log.actuator.command_limited(roll_rows,evaluation), [], "all")];
metrics.pitch_command_range = [ ...
    min(log.actuator.command_limited(pitch_rows,evaluation), [], "all"), ...
    max(log.actuator.command_limited(pitch_rows,evaluation), [], "all")];

metrics.average_power_W = mean(log.power.estimated_total_W(evaluation));
metrics.maximum_power_W = max(log.power.estimated_total_W(evaluation));

wrench_scale = max(abs(P.metrics.wrench_scale(:)), 1e-9);
normalized_residual = log.wrench.residual_B(:,evaluation)./wrench_scale;
metrics.wrench_residual_rms = ...
    sqrt(mean(log.wrench.residual_B(:,evaluation).^2,2));
metrics.normalized_wrench_residual_rms = ...
    sqrt(mean(normalized_residual.^2,2));

metrics.lost_control = numerical_divergence || ...
    any(position_error_norm > P.simulation.lost_control_position_error);

end


%% Quaternion and rotation utilities (no Aerospace Toolbox required)

function q = normalize_quaternion(q)

q = q(:);
n = norm(q);
if n < 1e-12 || ~isfinite(n)
    q = [1; 0; 0; 0];
else
    q = q/n;
end

end


function q_conjugate = quaternion_conjugate(q)

q = q(:);
q_conjugate = [q(1); -q(2:4)];

end


function q = quaternion_multiply(q1, q2)

q1 = q1(:);
q2 = q2(:);
s1 = q1(1);
v1 = q1(2:4);
s2 = q2(1);
v2 = q2(2:4);

q = [ ...
    s1*s2 - dot(v1,v2);
    s1*v2 + s2*v1 + cross(v1,v2)];

end


function R = quaternion_to_rotation(q)

q = normalize_quaternion(q);
w = q(1);
x = q(2);
y = q(3);
z = q(4);

R = [ ...
    1 - 2*(y^2 + z^2), 2*(x*y - w*z),     2*(x*z + w*y);
    2*(x*y + w*z),     1 - 2*(x^2 + z^2), 2*(y*z - w*x);
    2*(x*z - w*y),     2*(y*z + w*x),     1 - 2*(x^2 + y^2)];

end


function q = rotation_to_quaternion(R)

trace_R = trace(R);

if trace_R > 0
    s = 2*sqrt(trace_R + 1);
    q = [ ...
        0.25*s;
        (R(3,2) - R(2,3))/s;
        (R(1,3) - R(3,1))/s;
        (R(2,1) - R(1,2))/s];
elseif R(1,1) > R(2,2) && R(1,1) > R(3,3)
    s = 2*sqrt(max(1 + R(1,1) - R(2,2) - R(3,3), 0));
    q = [ ...
        (R(3,2) - R(2,3))/max(s,eps);
        0.25*s;
        (R(1,2) + R(2,1))/max(s,eps);
        (R(1,3) + R(3,1))/max(s,eps)];
elseif R(2,2) > R(3,3)
    s = 2*sqrt(max(1 + R(2,2) - R(1,1) - R(3,3), 0));
    q = [ ...
        (R(1,3) - R(3,1))/max(s,eps);
        (R(1,2) + R(2,1))/max(s,eps);
        0.25*s;
        (R(2,3) + R(3,2))/max(s,eps)];
else
    s = 2*sqrt(max(1 + R(3,3) - R(1,1) - R(2,2), 0));
    q = [ ...
        (R(2,1) - R(1,2))/max(s,eps);
        (R(1,3) + R(3,1))/max(s,eps);
        (R(2,3) + R(3,2))/max(s,eps);
        0.25*s];
end

q = normalize_quaternion(q);
if q(1) < 0
    q = -q;
end

end


function phi = rotation_log_vector(R)

cos_theta = 0.5*(trace(R) - 1);
cos_theta = min(max(cos_theta, -1), 1);
theta = acos(cos_theta);
skew_vector = 0.5*[ ...
    R(3,2) - R(2,3);
    R(1,3) - R(3,1);
    R(2,1) - R(1,2)];

if theta < 1e-7
    phi = skew_vector;
else
    phi = (theta/sin(theta))*skew_vector;
end

end


function limited = component_limit(value, absolute_limit)

absolute_limit = absolute_limit(:);
value = value(:);
if isscalar(absolute_limit)
    absolute_limit = repmat(absolute_limit, numel(value), 1);
end
limited = min(max(value, -absolute_limit), absolute_limit);

end
