function P = comet_parameters()
%COMET_PARAMETERS Parameters and experimental targets for the COMET model.
%
% The structure separates four kinds of information:
%   1) values reported in the COMET paper;
%   2) values derived from reported data;
%   3) transparent initial estimates required to run the simulation;
%   4) bounds used later for parameter calibration.
%
% Initial estimates must not be presented as measured COMET parameters.
% The calibration routine will update only the fields explicitly listed in
% P.calibration.

%% Metadata

P.meta.model_name = "COMET single-versus-dual swashplate model";
P.meta.version = "rewrite-v8-C2-entry-ramp-sweep";
P.meta.primary_controller = "paper_like_cascade";

P.meta.known_ambiguities = [
    "The paper reports z/D = 0.17 and z = 79 mm, while a nominal 16-inch propeller gives z/D = 0.194."
    "The inertia matrix, controller gains, derivative-filter bandwidths, actuator bandwidths and absolute swashplate effectiveness are not reported."
    "The paper calls its maximum absolute tracking error MAE; this model stores it as max_position_error."
    "Table III contains cyclic values rounded to +/-1 at 5 m/s although the text states that the dual configuration does not saturate."
];

%% Reference-frame conventions
%
% R maps body-frame vectors into the inertial frame:
%       vector_I = R * vector_B
%
% The inertial z-axis points upward. The nominal rotor thrust direction is
% +z_B, while gravity acts along -z_I.

P.frames.rotation_definition = "body_to_inertial";
P.frames.inertial_z_direction = "up";
P.frames.gravity_axis_I = [0; 0; -1];
P.frames.nominal_thrust_axis_B = [0; 0; 1];

%% Physical parameters

P.physical.g = 9.81;                         % [m/s^2] standard value
P.physical.mass = 1.340;                     % [kg] paper: tracking tests
P.physical.mass_tracking = 1.340;            % [kg] paper
P.physical.mass_endurance_standard = 1.250;  % [kg] paper: endurance test

P.physical.hover_thrust_total = ...
    P.physical.mass * P.physical.g;           % [N] derived
P.physical.hover_thrust_per_rotor = ...
    0.5 * P.physical.hover_thrust_total;      % [N] equal-thrust estimate

%% Geometry

P.geometry.height = 0.450;                    % [m] paper
P.geometry.folded_diameter = 0.066;           % [m] paper
P.geometry.body_radius = ...
    0.5 * P.geometry.folded_diameter;         % [m] simplified body model

P.geometry.propeller_diameter_nominal = ...
    16 * 0.0254;                              % [m] 16-inch propeller
P.geometry.projected_diameter = 0.465;        % [m] paper
P.geometry.effective_rotor_diameter = ...
    P.geometry.projected_diameter;            % [m] initial modelling choice
P.geometry.effective_rotor_radius = ...
    0.5 * P.geometry.effective_rotor_diameter;

P.geometry.rotor_separation = 0.079;          % [m] paper
P.geometry.z_over_D_reported = 0.17;          % [-] paper
P.geometry.z_over_nominal_propeller = ...
    P.geometry.rotor_separation / ...
    P.geometry.propeller_diameter_nominal;    % [-] derived consistency check

P.geometry.lower_blade_installation_angle = ...
    1.0 * pi / 180;                           % [rad] paper
P.geometry.flapping_phase_lag = pi / 4;       % [rad] paper: 45 deg

% Initial CoG assumption: midway between the two rotor hubs. These vectors
% must be updated if a better CoG estimate becomes available.
P.geometry.hub_position_up_B = ...
    [0; 0;  0.5 * P.geometry.rotor_separation];
P.geometry.hub_position_dw_B = ...
    [0; 0; -0.5 * P.geometry.rotor_separation];

%% Initial inertia estimate
%
% The paper does not report J. The estimate below uses a vertical solid
% cylinder for the central body and slender in-plane rods for the combined
% propeller mass. It is only an initialization for calibration/sensitivity.

P.physical.propeller_mass_assumed_total = 0.058; % [kg] paper value, assumed total

m_prop = P.physical.propeller_mass_assumed_total;
m_body = P.physical.mass - m_prop;
r_body = P.geometry.body_radius;
h_body = P.geometry.height;
r_rotor = P.geometry.effective_rotor_radius;

Jxx_body = (1/12) * m_body * (3*r_body^2 + h_body^2);
Jyy_body = Jxx_body;
Jzz_body = 0.5 * m_body * r_body^2;

Jxx_rotor = (1/6) * m_prop * r_rotor^2;
Jyy_rotor = Jxx_rotor;
Jzz_rotor = (1/3) * m_prop * r_rotor^2;

P.physical.J_initial = diag([ ...
    Jxx_body + Jxx_rotor, ...
    Jyy_body + Jyy_rotor, ...
    Jzz_body + Jzz_rotor]);                   % [kg m^2] estimated
P.physical.J = P.physical.J_initial;

%% Reported propulsion and hover data

P.propulsion.motor_model = "DT-5008";        % paper
P.propulsion.servo_model = "BlueArrow BA8-3";% paper
P.propulsion.propeller_model = "EOLO CN16x6";% paper

P.propulsion.hover_power_standard = 163.2;    % [W] at 1.250 kg, paper
P.propulsion.rotor_frequency_hover = 45.0;    % [Hz] paper vibration discussion
P.propulsion.rotor_speed_hover = ...
    2*pi*P.propulsion.rotor_frequency_hover;  % [rad/s] derived

% Initial drag-torque coefficient Q = k_Q*T, inferred from hover power.
% This assumes 45 Hz is the mechanical rotor frequency and equal loading.
T_hover_endurance = ...
    P.physical.mass_endurance_standard * P.physical.g;
P.aerodynamics.yaw_torque_per_thrust = ...
    P.propulsion.hover_power_standard / ...
    (P.propulsion.rotor_speed_hover * T_hover_endurance); % [m] derived estimate

%% Effective cyclic moment model
%
% Roll/pitch moment from rotor i will initially be modelled as
%       tau_i = c_i * T_i * delta_i
% where c_i has units of metres and delta_i is normalized in [-1,1].
%
% The absolute upper value is an INITIAL ESTIMATE. The lower/upper ratio and
% dual interaction factor are inferred from the reported maximum-torque
% improvements: dual/upper = 2.167 and dual/lower = 1.7251.

P.aerodynamics.dual_to_upper_torque_ratio = 2.167;   % paper
P.aerodynamics.dual_to_lower_torque_ratio = 1.7251;  % paper
P.aerodynamics.lower_to_upper_effectiveness = ...
    P.aerodynamics.dual_to_upper_torque_ratio / ...
    P.aerodynamics.dual_to_lower_torque_ratio;       % [-] derived

P.aerodynamics.dual_interaction_factor = ...
    P.aerodynamics.dual_to_upper_torque_ratio / ...
    (1 + P.aerodynamics.lower_to_upper_effectiveness); % [-] derived

P.aerodynamics.cyclic_arm_roll_up = 0.045;            % [m] INITIAL ESTIMATE
P.aerodynamics.cyclic_arm_pitch_up = 0.045;           % [m] INITIAL ESTIMATE

P.aerodynamics.cyclic_arm_roll_dw = ...
    P.aerodynamics.lower_to_upper_effectiveness * ...
    P.aerodynamics.cyclic_arm_roll_up;                 % [m] derived initial value
P.aerodynamics.cyclic_arm_pitch_dw = ...
    P.aerodynamics.lower_to_upper_effectiveness * ...
    P.aerodynamics.cyclic_arm_pitch_up;                % [m] derived initial value

% Global correction applied to every physical cyclic deflection before the
% thrust-vector and residual flapping moments are evaluated. A value below
% one reduces the command-to-wrench effectiveness while preserving the
% reported relative authority of the single and dual configurations.
P.aerodynamics.cyclic_effectiveness_scale = 0.780;       % [-] selected calibration value

%% Tip-path-plane mapping
%
% alpha_i = A_i * delta_ele_i
% beta_i  = B_i * delta_ail_i
%
% A_i and B_i are not reported. Five degrees at full normalized cyclic is a
% transparent INITIAL ESTIMATE used for force-vector visualization and the
% first simulation pass. These four values are calibration candidates.

tpp_angle_initial = 5.0 * pi / 180;                    % [rad] INITIAL ESTIMATE

P.aerodynamics.A_up = tpp_angle_initial;               % [rad/normalized input]
P.aerodynamics.B_up = tpp_angle_initial;
P.aerodynamics.A_dw = tpp_angle_initial;
P.aerodynamics.B_dw = tpp_angle_initial;

%% Aerodynamic drag
%
% The paper includes rotor/fuselage drag but does not identify its
% coefficients. These are INITIAL ESTIMATES and will be calibrated.

P.aerodynamics.linear_drag_B = ...
    diag([0.08, 0.08, 0.12]);                          % [N/(m/s)] estimate
P.aerodynamics.quadratic_drag_B = ...
    diag([0.0525, 0.0525, 0.00]);                     % [N/(m/s)^2] selected C2 value
P.aerodynamics.angular_damping_B = ...
    diag([0.002, 0.002, 0.001]);                       % [N m/(rad/s)] estimate

%% Actuator limits and dynamics
%
% Only the normalized cyclic interval is stated in the paper. Remaining
% actuator values are INITIAL ESTIMATES to be replaced by datasheet values
% or calibrated ranges.

P.actuators.cyclic_min = -1.0;                         % [-] paper
P.actuators.cyclic_max =  1.0;                         % [-] paper
P.actuators.cyclic_rate_max = 8.0;                     % [normalized/s] estimate
P.actuators.servo_time_constant = 0.020;               % [s] final rebalanced value

P.actuators.thrust_min = 0.0;                          % [N] physical constraint
P.actuators.thrust_max_per_rotor = 18.0;               % [N] estimate
P.actuators.thrust_rate_max = 100.0;                   % [N/s] estimate
P.actuators.motor_time_constant = 0.080;               % [s] estimate

P.actuators.apply_position_limits = true;
P.actuators.apply_rate_limits = true;
P.actuators.saturation_tolerance = 1e-6;

%% Configuration masks
%
% Cyclic order:
%   [ail_up; ele_up; ail_dw; ele_dw]
%
% In the paper's single-swashplate flight test, the upper cyclic mechanism
% is disabled and the lower one remains active.

P.configuration.single.cyclic_mask = [0; 0; 1; 1];
P.configuration.dual.cyclic_mask = [1; 1; 1; 1];
P.configuration.single.name = "single-lower";
P.configuration.dual.name = "dual";

%% Mixer and control-allocation parameters
%
% A common roll/pitch command is sent to every active swashplate. The mixer
% yaw scale below is common to both configurations. The nominal dual roll
% and pitch scales are retained only for metric normalization; the actual
% allocation uses a configuration-dependent local effectiveness matrix.

T_h = P.physical.hover_thrust_per_rotor;

P.mixer.roll_moment_per_command = ...
    P.aerodynamics.dual_interaction_factor * T_h * ...
    (P.aerodynamics.cyclic_arm_roll_up + ...
     P.aerodynamics.cyclic_arm_roll_dw);                % [N m]

P.mixer.pitch_moment_per_command = ...
    P.aerodynamics.dual_interaction_factor * T_h * ...
    (P.aerodynamics.cyclic_arm_pitch_up + ...
     P.aerodynamics.cyclic_arm_pitch_dw);               % [N m]

P.mixer.yaw_moment_per_thrust_difference = ...
    P.aerodynamics.yaw_torque_per_thrust;                % [m]

P.mixer.independent_swashplate_commands = false;

% At every controller sample, simulate_COMET linearizes the complete rotor
% wrench model around the requested rotor thrusts. This includes both the
% residual flapping moments and the geometric r_hub x F terms. Consequently,
% the single configuration must request a larger normalized cyclic command
% to produce the same desired roll/pitch moment.
P.allocation.configuration_dependent_cyclic_map = true;
P.allocation.cyclic_linearization_step = 1e-4;           % normalized command
P.allocation.minimum_reciprocal_condition = 1e-8;
P.allocation.use_pseudoinverse_if_ill_conditioned = true;

%% Paper-like cascaded controller: initial gains
%
% Numerical gains are not reported in the paper. These values are initial
% guesses informed by the previous stable simplified model. They must be
% tuned once and then kept identical for single and dual configurations.

P.controller.type = "paper_like_cascade";

% Position loop: position error -> desired velocity correction.
P.controller.position.Kp = diag([1.0, 1.0, 1.5]);
P.controller.position.velocity_command_limit = [6; 6; 3]; % [m/s]

% Velocity PID plus acceleration feedforward.
P.controller.velocity.Kp = diag([3.0, 3.0, 4.0]);
P.controller.velocity.Ki = diag([0.10, 0.10, 0.15]);
P.controller.velocity.Kd = diag([0.15, 0.15, 0.20]);
P.controller.velocity.derivative_filter_time_constant = 0.050; % [s] estimate
P.controller.velocity.integral_limit = [2; 2; 2];
P.controller.velocity.acceleration_limit = [12; 12; 8]; % [m/s^2]

% Quaternion attitude outer loop: attitude error -> desired angular rate.
P.controller.attitude.Kp = diag([4.0, 4.0, 2.0]);
P.controller.attitude.omega_command_limit = [8; 8; 5];  % [rad/s]

% Angular-rate PID: angular-rate error -> desired body moment.
P.controller.rate.Kp = diag([0.45, 0.45, 0.25]);
P.controller.rate.Ki = diag([0.03, 0.03, 0.02]);
P.controller.rate.Kd = diag([0.015, 0.015, 0.010]);
P.controller.rate.derivative_filter_time_constant = 0.030; % [s] estimate
P.controller.rate.integral_limit = [1; 1; 1];

P.controller.anti_windup_gain = 1.0;
P.controller.same_gains_for_all_configurations = true;

%% Optional geometric-controller extension
%
% These reproduce the previous project's initial gains and are retained
% only for the later extension. They are not used in the paper baseline.

P.controller.geometric.Kp_position = diag([4.0, 4.0, 6.0]);
P.controller.geometric.Kv = diag([3.0, 3.0, 4.0]);
P.controller.geometric.KR = diag([2.0, 2.0, 1.0]);
P.controller.geometric.Komega = diag([0.4, 0.4, 0.3]);

%% Trajectory geometry

P.trajectory.hover_height = 1.0;                      % [m]
P.trajectory.circle_radius = 2.5;                     % [m] paper: 5 m diameter
P.trajectory.eight_A = 5.0;                           % [m] x range +/-5 m
P.trajectory.eight_B = 2.5;                           % [m] y range +/-2.5 m
P.trajectory.default_yaw = 0.0;                       % [rad]
P.trajectory.default_cycles = 2;
P.trajectory.start_from_rest = true;
P.trajectory.use_acceleration_limited_entry = true;
% Final duration selected from the entry-ramp sensitivity study.
P.trajectory.entry_time_override = 1.30;               % [s]

%% Simulation settings

% The digital controller is always evaluated at 200 Hz. RK4 may use one or
% more internal substeps while the controller output is held constant. This
% separates a numerical-convergence test from a change in controller rate.
P.simulation.controller_sample_time = 0.005;           % [s] fixed 200 Hz
P.simulation.integration_substeps = 1;                 % [-] nominal RK4
P.simulation.dt = P.simulation.controller_sample_time; % [s] legacy alias
P.simulation.dt_convergence_check = ...
    P.simulation.controller_sample_time/2;             % [s] two-substep RK4
P.simulation.integrator = "RK4";
P.simulation.normalize_quaternion_each_step = true;

P.simulation.max_position_norm = 30.0;                % [m]
P.simulation.max_velocity_norm = 20.0;                % [m/s]
P.simulation.max_angular_rate = 80.0;                 % [rad/s]
P.simulation.lost_control_position_error = 2.0;       % [m] analysis threshold
P.simulation.stop_on_numerical_divergence = true;
P.simulation.stop_on_lost_control = false;

%% Measurement and feedback delay

% Physical delay used by the feedback controller. simulate_COMET applies
% fractional-sample delays by interpolation, so 2.5 ms is not rounded to
% the 5 ms controller sample time.
P.feedback.state_delay = 0.0025;                       % [s] final physical value

% Timestamp lag used only when run_COMET_comparison evaluates the reported
% position-tracking metrics. It is not used by the feedback controller and
% therefore cannot change dynamics, power or actuator saturation.
P.measurement.position_timestamp_lag = 0.060;          % [s] selected diagnostic value

%% Disturbance defaults

P.disturbance.enabled = false;
P.disturbance.type = "none";

P.disturbance.gust.force_I = [2.0; 0; 0];             % [N] initial test value
P.disturbance.gust.start_time = 5.0;                   % [s]
P.disturbance.gust.duration = 1.0;                     % [s]
P.disturbance.gust.profile = "one_cosine";

P.disturbance.torque.value_B = [0.15; 0; 0];          % [N m]
P.disturbance.torque.start_time = 5.0;                 % [s]
P.disturbance.torque.duration = 0.25;                  % [s]

P.disturbance.random_seed = 1968652;
P.disturbance.parameter_uncertainty_levels = [0.10, 0.20];

%% Metrics

tau_roll_scale = P.mixer.roll_moment_per_command;
tau_pitch_scale = P.mixer.pitch_moment_per_command;
tau_yaw_scale = ...
    P.aerodynamics.yaw_torque_per_thrust * ...
    P.physical.hover_thrust_total;

P.metrics.wrench_scale = [ ...
    P.physical.hover_thrust_total;
    P.physical.hover_thrust_total;
    P.physical.hover_thrust_total;
    tau_roll_scale;
    tau_pitch_scale;
    tau_yaw_scale];

% Include the initial transient in RMSE, maximum error and saturation metrics.
P.metrics.evaluate_after_entry = false;
P.metrics.report_raw_and_normalized_wrench_residual = true;

%% Calibration variables and bounds
%
% Only these values may be changed automatically in the first calibration
% stage. Controller gains will be tuned separately to avoid compensating for
% an incorrect physical model.

P.calibration.names = [
    "cyclic_arm_roll_up"
    "cyclic_arm_pitch_up"
    "lower_to_upper_effectiveness"
    "dual_interaction_factor"
    "servo_time_constant"
    "linear_drag_xy"
];

P.calibration.initial = [
    P.aerodynamics.cyclic_arm_roll_up
    P.aerodynamics.cyclic_arm_pitch_up
    P.aerodynamics.lower_to_upper_effectiveness
    P.aerodynamics.dual_interaction_factor
    P.actuators.servo_time_constant
    P.aerodynamics.linear_drag_B(1,1)
];

P.calibration.lower_bound = [
    0.010
    0.010
    1.000
    0.800
    0.010
    0.000
];

P.calibration.upper_bound = [
    0.150
    0.150
    1.600
    1.100
    0.150
    0.500
];

%% Experimental targets from the paper

P.experimental.tracking = paper_tracking_cases();

P.experimental.hover.energy_Wh = [66.6, 111.0, 111.0];
P.experimental.hover.takeoff_mass_kg = [1.250, 1.460, 1.940];
P.experimental.hover.power_W = [163.2, 204.0, 294.0];
P.experimental.hover.time_min = [24.3, 35.5, 22.0];
P.experimental.hover.efficiency_g_per_W = [7.66, 7.15, 6.00];

P.experimental.force_sensor.max_force_error = 0.04;   % [N] paper
P.experimental.force_sensor.max_torque_error = 0.009; % [N m] paper
P.experimental.bench_test.throttle = 0.40;            % [-] paper

%% Parameter provenance labels

P.status.mass = "reported";
P.status.geometry = "reported_with_D_ambiguity";
P.status.inertia = "estimated";
P.status.torque_ratios = "reported";
P.status.absolute_cyclic_effectiveness = "initial_estimate_to_calibrate";
P.status.cyclic_effectiveness_scale = "selected_calibration_value";
P.status.tpp_mapping = "initial_estimate_to_calibrate";
P.status.drag = "selected_calibration_value";
P.status.actuator_dynamics = "selected_rebalanced_value";
P.status.feedback_delay = "selected_rebalanced_value";
P.status.entry_ramp = "selected_calibration_value";
P.status.measurement_timestamp_lag = "selected_diagnostic_value";
P.status.controller_gains = "initial_estimate_to_tune_once";
P.status.derivative_filters = "initial_estimate_to_tune_once";

%% Internal consistency checks

validate_parameter_structure(P);

end

%% Local helper: paper trajectory-tracking targets

function cases = paper_tracking_cases()

cases = repmat(empty_case(), 1, 5);

cases(1) = make_case("circle_3", "circle", 3.0, 3.71);
cases(1).single = make_target(true, true, 0.2908, 0.6161, ...
    [-0.8013, 0.5770], [-0.6767, 0.9912], 187.98, 223.22);
cases(1).dual = make_target(true, true, 0.1712, 0.3417, ...
    [-0.3363, 0.5214], [-0.3984, 0.6028], 170.57, 183.75);

cases(2) = make_case("circle_4", "circle", 4.0, 6.64);
cases(2).single = make_target(true, true, 0.3924, 0.8696, ...
    [-1.0000, 0.9893], [-0.8585, 1.0000], 201.86, 262.49);
cases(2).dual = make_target(true, true, 0.3139, 0.6201, ...
    [-0.6473, 0.5900], [-0.5806, 0.6938], 198.39, 257.99);

cases(3) = make_case("eight_3p5", "eight", 3.5, 5.00);
cases(3).single = make_target(true, true, 0.3221, 0.7833, ...
    [-1.0000, 0.8318], [-0.9308, 0.4953], 194.78, 246.48);
cases(3).dual = make_target(true, true, 0.2534, 0.5165, ...
    [-0.7554, 0.5490], [-0.6379, 0.4239], 186.93, 215.85);

cases(4) = make_case("eight_4", "eight", 4.0, 6.45);
cases(4).single = make_target(true, false, NaN, NaN, ...
    [NaN, NaN], [NaN, NaN], NaN, NaN);
cases(4).single.note = "Lost control at the first turn due to actuator saturation.";
cases(4).dual = make_target(true, true, 0.2876, 0.5650, ...
    [-0.8821, 0.3635], [-0.4064, 0.7926], 193.13, 236.77);

cases(5) = make_case("eight_5", "eight", 5.0, 6.79);
cases(5).single = make_target(false, false, NaN, NaN, ...
    [NaN, NaN], [NaN, NaN], NaN, NaN);
cases(5).single.note = "Single-swashplate result not reported.";
cases(5).dual = make_target(true, true, 0.3052, 0.6526, ...
    [-1.0000, 0.4353], [-0.5734, 1.0000], 197.08, 237.03);

end


function test = empty_case()

test.id = "";
test.trajectory = "";
test.vmax = NaN;
test.amax = NaN;
test.single = empty_target();
test.dual = empty_target();

end


function test = make_case(id, trajectory, vmax, amax)

test = empty_case();
test.id = string(id);
test.trajectory = string(trajectory);
test.vmax = vmax;
test.amax = amax;

end


function target = empty_target()

target.available = false;
target.completed = false;
target.rmse_position = NaN;
target.max_position_error = NaN;
target.roll_command_range = [NaN, NaN];
target.pitch_command_range = [NaN, NaN];
target.average_power = NaN;
target.maximum_power = NaN;
target.note = "";

end


function target = make_target(available, completed, rmse, max_error, ...
    roll_range, pitch_range, average_power, maximum_power)

target = empty_target();
target.available = logical(available);
target.completed = logical(completed);
target.rmse_position = rmse;
target.max_position_error = max_error;
target.roll_command_range = roll_range;
target.pitch_command_range = pitch_range;
target.average_power = average_power;
target.maximum_power = maximum_power;

end


function validate_parameter_structure(P)

assert(P.physical.mass > 0, "Vehicle mass must be positive.");
assert(all(diag(P.physical.J) > 0), ...
    "The diagonal inertia entries must be positive.");
assert(P.geometry.rotor_separation > 0, ...
    "Rotor separation must be positive.");
assert(P.actuators.cyclic_min < P.actuators.cyclic_max, ...
    "Invalid cyclic limits.");
assert(P.actuators.thrust_max_per_rotor > ...
       P.physical.hover_thrust_per_rotor, ...
    "Estimated maximum rotor thrust must exceed hover thrust.");
assert(P.simulation.controller_sample_time > 0, ...
    "Controller sample time must be positive.");
assert(P.feedback.state_delay >= 0 && ...
       isfinite(P.feedback.state_delay), ...
    "The feedback-state delay must be finite and nonnegative.");
assert(P.measurement.position_timestamp_lag >= 0 && ...
       isfinite(P.measurement.position_timestamp_lag), ...
    "The measurement timestamp lag must be finite and nonnegative.");
assert(isnan(P.trajectory.entry_time_override) || ...
       (isscalar(P.trajectory.entry_time_override) && ...
        isfinite(P.trajectory.entry_time_override) && ...
        P.trajectory.entry_time_override > 0), ...
    "The entry-time override must be NaN or a positive finite scalar.");
assert(P.simulation.integration_substeps >= 1 && ...
       P.simulation.integration_substeps == ...
       round(P.simulation.integration_substeps), ...
    "The RK4 integration-substep count must be a positive integer.");
assert(P.controller.velocity.derivative_filter_time_constant > 0, ...
    "The velocity derivative-filter time constant must be positive.");
assert(P.controller.rate.derivative_filter_time_constant > 0, ...
    "The rate derivative-filter time constant must be positive.");
assert(P.allocation.cyclic_linearization_step > 0 && ...
       P.allocation.cyclic_linearization_step < 0.1, ...
    "The cyclic-allocation linearization step is invalid.");
assert(P.allocation.minimum_reciprocal_condition > 0, ...
    "The minimum reciprocal condition must be positive.");
assert(P.aerodynamics.cyclic_effectiveness_scale > 0 && ...
       P.aerodynamics.cyclic_effectiveness_scale <= 1, ...
    "The cyclic-effectiveness scale must lie in the interval (0,1].");
assert(all(P.calibration.lower_bound < P.calibration.upper_bound), ...
    "Every calibration lower bound must be smaller than its upper bound.");
assert(numel(P.calibration.names) == numel(P.calibration.initial), ...
    "Calibration names and initial values have inconsistent sizes.");

end
