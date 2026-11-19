COMET — INDEPENDENT HYBRID DIFFERENTIAL-FLATNESS EXPERIMENT

Keep the original three files unchanged in your MATLAB folder:
  comet_parameters.m
  simulate_COMET.m
  run_COMET_comparison.m

Add ONLY:
  simulate_COMET_flatness.m
  run_COMET_flatness.m

Run the original study with run_COMET_comparison.
Run the extended independent study with run_COMET_flatness.

The independent runner performs ten baseline runs and ten hybrid-flatness runs.
It retains the original statistics and plots, then adds comparisons and
nominal thrust/moment diagnostics. The original runner is NEVER modified.

Architecture:
1. Existing analytic circle/eight position and yaw references.
2. NEW analytic derivatives through jerk and snap (including the entry ramp).
3. Reduced quadrotor-flatness inverse dynamics: nominal thrust T,
   nominal rotation R, desired angular rate Omega, angular acceleration
   OmegaDot, nominal body moment J*OmegaDot + cross(Omega,J*Omega).
4. EXISTING cascaded position/velocity feedback computes corrected force
   and desired attitude. Existing attitude and angular-rate feedback remain.
5. Nominal angular rate is transported into actual body coordinates;
   feedforward angular acceleration includes the required rotating-frame
   transport term before it is added to the rate-feedback torque.
6. EXISTING COMET nonlinear 19-state plant, dual/single masking, rotor
   force/moment map, control allocation, actuator dynamics, limits, RK4,
   delayed feedback and metric definitions are retained.

IMPORTANT SCOPE:
- The flatness reconstruction is for a near-hover, quadrotor-like REDUCED
  model. It does NOT reconstruct COMET's six physical actuator states or
  servo/motor commands from flat outputs.
- The physical COMET control mixer determines feasible actuator commands.
- Thrust nominal T is computed and logged. The actual total thrust demand
  continues to include feedback corrections; it is NOT blindly replaced
  with open-loop flatness thrust.
- The nominal moment is computed and logged. Its angular-acceleration
  component, expressed in current body coordinates with transport, is
  combined with the existing angular-rate PID and gyroscopic compensation.
- The original trajectory entry is a quintic VELOCITY ramp; it is not
  globally snap-continuous at the entry/constant-speed junction. This
  can cause jumps in the nominal angular-acceleration feedforward there.
- Not a proof that the full 19-state COMET is differentially flat.
- Unlike the original paper's explicitly documented planner outputs,
  use of fourth derivatives and moment feedforward here is an extension.

CHECKS AFTER RUN:
- Confirm the nominal thrust and moment diagnostics are finite.
- Inspect discontinuities near entry_time=1.30 s.
- Inspect moment and cyclic saturation and RMSE for both architectures.
- Baseline results should remain numerically the same as the original run.
- MATLAB was NOT available in the authoring environment. Integration and
  numerical behavior have NOT been validated by executing MATLAB.
