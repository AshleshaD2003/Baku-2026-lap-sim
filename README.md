# 2D Point-Mass Lap Time Simulator — Baku City Circuit (2026 F1 Spec)

A quasi-steady-state 2D point-mass lap time simulator written in Python, configured for the **Baku City Circuit** using the **2026 Formula 1 Technical Regulations**. 

The model parses vector circuit geometry from an SVG file, computes true Euclidean waypoint spacing and local path curvature, and solves for the maximum velocity profile across forward-acceleration, backward braking, and cornering limits.

**SVG Path Parsing & Re-scaling:** Extracts cubic Bézier curves from circuit outlines, scales the path to Baku's official 6.003 km length, and resamples waypoints.Euclidean Kinematics Fix: Calculates individual variable step sizes between discrete track coordinates to eliminate arc-length distortion during forward/backward numerical integration.

**Wide-Baseline Curvature Fitting:** Implements a 3-point circumradius algorithm across a wide spatial baseline to filter SVG digitisation noise while capturing accurate corner radii.

**2026 Active Aerodynamics:** Dynamic switching between X-Mode (High Downforce) in cornering zones and Z-Mode (Low Drag) on straights.

**Hybrid Powertrain & Energy Management:** Models a ~400 kW Internal Combustion Engine (ICE) alongside a ~350 kW MGU-K with a per-lap kinetic energy recovery limit (8.5MJ) and friction-circle traction constraints.

**Environmental & Grip Sensitivity:** Built-in track condition modelling adjusting air density and tyre friction for ambient temperatures and weather conditions (Dry, Damp, Wet).

**Model Assumptions**

1. Point-Mass Assumption: Vehicle rotational dynamics, lateral weight transfer, roll center, and suspension kinematics are simplified to a point-mass friction circle.

2. Auto-Numbered Geometry: Track corners are evaluated sequentially based on path sampling rather than official turn designations.

3. Manual Override Mode (MOM): 2026 MOM boost mode is turned off by default for single-car baseline qualifying simulation.
