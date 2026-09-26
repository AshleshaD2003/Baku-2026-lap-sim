# 2D Point-Mass Lap Time Simulator — Baku City Circuit (2026 F1 Spec)

A quasi-steady-state 2D point-mass lap time simulator written in Python, configured for the **Baku City Circuit** using the **2026 Formula 1 Technical Regulations**. 

The model parses vector circuit geometry from an SVG file, computes true Euclidean waypoint spacing and local path curvature, and solves for the maximum velocity profile across forward-acceleration, backward-braking, and cornering limits.