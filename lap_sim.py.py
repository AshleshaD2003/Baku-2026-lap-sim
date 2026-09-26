#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().system('pip install numpy matplotlib pandas')
"""
Baku City Circuit — Point-Mass Lap Time Simulation (2026-spec F1 car)
======================================================================

Pipeline:
  1. Parse the traced circuit SVG into a dense (x, y) polyline.
  2. Scale it to the real centreline length (6.003 km) and resample at a
     nominal spacing.
  3. Compute actual point-by-point Euclidean step distances (real_ds) to eliminate
     arc-length distortion errors in kinematics and time integration.
  4. Compute local curvature with a wide-baseline 3-point circle fit.
  5. Simulate with a friction-circle-limited point mass car using 2026 technical
     regulations (active aero, energy-limited MGU-K, variable track conditions).
"""

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# 1. SVG PATH PARSING & GEOMETRY RECONSTRUCTION
# ---------------------------------------------------------------------------

def parse_path_to_polyline(svg_path, samples_per_curve=25):
    """Flatten an SVG path 'd' string (m/c/s/z, relative commands) into a
    dense (x, y) polyline by sampling every cubic Bezier segment."""
    token_re = re.compile(r'[MmCcSsLlZz]|-?\d*\.?\d+(?:[eE][-+]?\d+)?')
    tokens = token_re.findall(svg_path)

    i = 0
    cur = np.array([0.0, 0.0])
    last_c2 = None
    last_cmd = None
    points = [cur.copy()]

    def read_nums(n):
        nonlocal i
        vals = [float(tokens[i + k]) for k in range(n)]
        i += n
        return vals

    def cubic_bezier(p0, c1, c2, p3, n):
        pts = []
        for t in np.linspace(0, 1, n)[1:]:
            pts.append(((1 - t) ** 3) * p0 + 3 * ((1 - t) ** 2) * t * c1
                       + 3 * (1 - t) * (t ** 2) * c2 + (t ** 3) * p3)
        return pts

    while i < len(tokens):
        tok = tokens[i]
        if tok in 'MmCcSsLlZz':
            cmd = tok
            i += 1
        else:
            cmd = last_cmd

        if cmd == 'm':
            dx, dy = read_nums(2)
            cur = cur + np.array([dx, dy])
            points.append(cur.copy())
            last_c2 = None
            last_cmd = 'l'
        elif cmd == 'c':
            dx1, dy1, dx2, dy2, dx, dy = read_nums(6)
            c1, c2 = cur + [dx1, dy1], cur + [dx2, dy2]
            p3 = cur + np.array([dx, dy])
            points.extend(cubic_bezier(cur, c1, c2, p3, samples_per_curve))
            cur, last_c2, last_cmd = p3, c2, 'c'
        elif cmd == 's':
            dx2, dy2, dx, dy = read_nums(4)
            c1 = cur + (cur - last_c2) if last_c2 is not None else cur.copy()
            c2 = cur + np.array([dx2, dy2])
            p3 = cur + np.array([dx, dy])
            points.extend(cubic_bezier(cur, c1, c2, p3, samples_per_curve))
            cur, last_c2, last_cmd = p3, c2, 's'
        elif cmd == 'l':
            dx, dy = read_nums(2)
            cur = cur + np.array([dx, dy])
            points.append(cur.copy())
            last_c2, last_cmd = None, 'l'
        elif cmd in 'zZ':
            points.append(points[0])
            last_cmd = None
        else:
            raise ValueError(f"Unhandled path command: {cmd}")

    pts = np.array(points)
    return pts[:, 0], pts[:, 1]


def load_track_from_svg(svg_file, real_length_m, ds_nominal=2.0):
    """Parse, scale to real length, and compute true per-step Euclidean distances."""
    with open(svg_file) as f:
        content = f.read()
    d = re.search(r']*\sd="([^"]+)"', content).group(1)
    x_raw, y_raw = parse_path_to_polyline(d)

    seg = np.sqrt(np.diff(x_raw) ** 2 + np.diff(y_raw) ** 2)
    cum = np.concatenate([[0], np.cumsum(seg)])
    scale = real_length_m / cum[-1]

    s_nominal = np.arange(0, real_length_m, ds_nominal)
    x = np.interp(s_nominal, cum * scale, x_raw * scale)
    y = np.interp(s_nominal, cum * scale, y_raw * scale)

    # Compute actual step-by-step Euclidean distances (real_ds)
    dx = np.diff(x)
    dy = np.diff(y)
    real_ds = np.sqrt(dx**2 + dy**2)
    s_actual = np.insert(np.cumsum(real_ds), 0, 0.0)

    return s_actual, x, y, real_ds


def curvature_wide_baseline(x, y, stride=8):
    """3-point circumradius curvature with a wide baseline to filter digitization noise."""
    n = len(x)
    curvature = np.zeros(n)
    for i in range(stride, n - stride):
        A = np.array([x[i - stride], y[i - stride]])
        B = np.array([x[i], y[i]])
        C = np.array([x[i + stride], y[i + stride]])
        ab, bc, ca = np.linalg.norm(B - A), np.linalg.norm(C - B), np.linalg.norm(A - C)
        cross = abs((B[0]-A[0])*(C[1]-B[1]) - (B[1]-A[1])*(C[0]-B[0]))
        curvature[i] = 0.0 if cross < 1e-9 else 1.0 / ((ab * bc * ca) / (2.0 * cross))
    curvature[:stride] = curvature[stride]
    curvature[-stride:] = curvature[-stride - 1]
    return curvature


# ---------------------------------------------------------------------------
# 2. TRACK CONDITIONS
# ---------------------------------------------------------------------------

@dataclass
class TrackConditions:
    air_temp_c: float = 26.0      
    track_temp_c: float = 42.0    
    pressure_pa: float = 101000.0 
    weather: str = "dry"          
    grip_level: float = 1.0       

    def air_density(self):
        return self.pressure_pa / (287.05 * (self.air_temp_c + 273.15))

    def mu_multiplier(self):
        weather_factor = {"dry": 1.00, "damp": 0.75, "wet": 0.55}[self.weather]
        temp_penalty = max(0.6, 1.0 - 0.0006 * (self.track_temp_c - 45.0) ** 2)
        return weather_factor * self.grip_level * temp_penalty


# ---------------------------------------------------------------------------
# 3. VEHICLE MODEL (2026 TECHNICAL REGULATIONS)
# ---------------------------------------------------------------------------

class F1Vehicle2026:
    def __init__(self, conditions: TrackConditions,
                 mass=768.0,                 
                 area=1.6,
                 mu_base=1.7,
                 ice_power=400_000,          
                 mgu_k_power=350_000,        
                 mgu_k_energy_per_lap=8.5e6, 
                 Cl_corner=3.6, Cd_corner=1.25,   
                 Cl_straight=1.8, Cd_straight=0.70):  
        self.mass = mass
        self.g = 9.81
        self.area = area
        self.mu_base = mu_base
        self.ice_power = ice_power
        self.mgu_k_power = mgu_k_power
        self.max_power = ice_power + mgu_k_power
        self.mgu_k_capacity = mgu_k_energy_per_lap
        self.Cl_corner, self.Cd_corner = Cl_corner, Cd_corner
        self.Cl_straight, self.Cd_straight = Cl_straight, Cd_straight
        self.conditions = conditions
        self.rho = conditions.air_density()
        self.mu = mu_base * conditions.mu_multiplier()
        self.battery_j = mgu_k_energy_per_lap  

    def aero_mode(self, radius):
        if radius > 400:
            return self.Cl_straight, self.Cd_straight
        return self.Cl_corner, self.Cd_corner

    def _forces(self, v, radius):
        v = max(v, 1.0)
        Cl, Cd = self.aero_mode(radius)
        downforce = 0.5 * self.rho * Cl * self.area * v**2
        drag = 0.5 * self.rho * Cd * self.area * v**2
        f_normal = self.mass * self.g + downforce
        f_max = self.mu * f_normal
        f_lat = 0.0 if radius > 2000 else min((self.mass * v**2) / radius, f_max)
        f_avail = np.sqrt(max(f_max**2 - f_lat**2, 0.0))
        return drag, f_avail

    def corner_speed_limit(self, radius):
        Cl, _ = self.aero_mode(radius)
        if radius > 2000:
            return 350.0 / 3.6
        num = self.mu * self.mass * self.g
        den = (self.mass / radius) - (0.5 * self.rho * Cl * self.area * self.mu)
        if den <= 0:
            return 350.0 / 3.6
        return np.sqrt(num / den)

    def max_acceleration(self, v, radius, ds_i, use_boost=False):
        drag, f_avail = self._forces(v, radius)
        v_eff = max(v, 1.0)

        available_power = self.ice_power
        if self.battery_j > 0:
            available_power += self.mgu_k_power
        if use_boost and self.battery_j > 0 and v < 337 / 3.6:
            available_power += self.mgu_k_power

        f_power = available_power / v_eff
        f_drive = min(f_avail, f_power)
        a_net = max((f_drive - drag) / self.mass, 0.0)

        if self.battery_j > 0 and f_drive > 0:
            elec_fraction = min(1.0, max(0.0, (f_drive * v_eff - self.ice_power) / max(self.mgu_k_power, 1)))
            self.battery_j = max(0.0, self.battery_j - elec_fraction * self.mgu_k_power * (ds_i / v_eff))
        return a_net

    def max_deceleration(self, v, radius, ds_i, regen_efficiency=0.9):
        drag, f_avail = self._forces(v, radius)
        regen_j = regen_efficiency * f_avail * ds_i
        self.battery_j = min(self.mgu_k_capacity, self.battery_j + regen_j * 0.05)
        return (f_avail + drag) / self.mass


# ---------------------------------------------------------------------------
# 4. LAP SIMULATION (Variable Step-Size Integration)
# ---------------------------------------------------------------------------

def find_straights(radii, s, threshold=400, min_len=100):
    is_straight = radii > threshold
    segs, i, n = [], 0, len(is_straight)
    while i < n:
        if is_straight[i]:
            j = i
            while j < n and is_straight[j]:
                j += 1
            length = s[j - 1] - s[i]
            if length >= min_len:
                segs.append((s[i], length))
            i = j
        else:
            i += 1
    return sorted(segs, key=lambda t: -t[1])


def run_lap(car, s, radii, real_ds):
    n = len(s)
    v_apex = np.array([car.corner_speed_limit(r) for r in radii])

    # Forward Pass: Accelerate step-by-step using individual ds_i
    v_forward = np.copy(v_apex)
    for i in range(1, n):
        ds_i = real_ds[i - 1]
        a = car.max_acceleration(v_forward[i - 1], radii[i - 1], ds_i)
        v_forward[i] = min(v_apex[i], np.sqrt(v_forward[i - 1]**2 + 2 * a * ds_i))

    # Backward Pass: Decelerate step-by-step using individual ds_i
    v_sim = np.copy(v_forward)
    for i in range(n - 2, -1, -1):
        ds_i = real_ds[i]
        a = car.max_deceleration(v_sim[i + 1], radii[i + 1], ds_i)
        v_sim[i] = min(v_forward[i], np.sqrt(v_sim[i + 1]**2 + 2 * a * ds_i))

    # Lap Time Integration: Match element length with v_sim
    ds_full = np.append(real_ds, real_ds[-1])
    lap_time = np.sum(ds_full / v_sim)
    return lap_time, v_sim, v_apex


# ---------------------------------------------------------------------------
# 5. EXECUTION & VISUALIZATION
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    REAL_LENGTH_M = 6003.0   # Baku official centreline length
    DS_NOMINAL = 2.0

    s_actual, x, y, real_ds = load_track_from_svg("RaceCircuitGillesBaku.svg", REAL_LENGTH_M, ds_nominal=DS_NOMINAL)
    curv = curvature_wide_baseline(x, y, stride=8)   
    radii = np.where(curv > 1e-6, 1.0 / curv, np.inf)

    straights = find_straights(radii, s_actual)
    print("Longest straight sections detected (start_m, length_m):")
    for st, ln in straights[:5]:
        print(f"  start={st:6.0f} m   length={ln:6.0f} m")

    conditions = TrackConditions(weather="dry", grip_level=1.0, air_temp_c=26.0, track_temp_c=42.0)
    car = F1Vehicle2026(conditions)

    lap_time, v_sim, v_apex = run_lap(car, s_actual, radii, real_ds)
    print(f"\n--- 2026-spec car, dry, rubbered track ---")
    print(f"Actual Track Length: {s_actual[-1]:.2f} m")
    print(f"Lap time: {lap_time:.3f} s   Top speed: {v_sim.max()*3.6:.1f} km/h   Slowest corner: {v_sim.min()*3.6:.1f} km/h\n")

    # Quick Sensitivity Check: Dry vs Wet
    for weather, grip in [("dry", 1.0), ("damp", 0.9), ("wet", 0.8)]:
        c2 = TrackConditions(weather=weather, grip_level=grip)
        car2 = F1Vehicle2026(c2)
        t2, v2, _ = run_lap(car2, s_actual, radii, real_ds)
        print(f"  {weather:>5s} track (grip_level={grip}): lap time = {t2:.2f} s ({(t2-lap_time)/lap_time*100:+.1f}% vs dry)")

    # --- Telemetry Plots ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc = axes[0].scatter(x, y, c=v_sim*3.6, cmap='RdYlGn', s=4)
    axes[0].set_aspect('equal')
    axes[0].set_title("Baku City Circuit — Speed by Position")
    plt.colorbar(sc, ax=axes[0], label="km/h")

    axes[1].plot(s_actual, v_sim*3.6, color='#00d2be', lw=1.5, label='Simulated speed')
    axes[1].plot(s_actual, np.minimum(v_apex*3.6, 400), color='#ff1801', lw=0.8, alpha=0.5, label='Apex limit (clipped)')
    axes[1].set_xlabel("Distance (m)")
    axes[1].set_ylabel("Speed (km/h)")
    axes[1].set_title(f"Speed Trace — Lap Time {lap_time:.2f}s (Dry)")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig("baku_telemetry.png", dpi=150)
    print("\nSaved plot to baku_telemetry.png")


# In[ ]:




