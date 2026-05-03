# Theory Note

## Interface Temperature Feedback

The control variable is the interface temperature measured on the face zone identified by `ID_INTERFACE`. The UDF samples only faces whose radial coordinate lies within:

$$
R_{\min} \le r \le R_{\max}
$$

with:

$$
r = \sqrt{x^2 + z^2}
$$

This restricts feedback to the intended annular region of the crystal/melt interface.

## Area-Weighted Temperature Measurement

For all selected interface faces, the UDF forms an area-weighted average temperature:

$$
T_{\mathrm{avg}} = \frac{\sum_i T_i A_i}{\sum_i A_i}
$$

where \(A_i\) is the face area magnitude and \(T_i\) is the Fluent face temperature. This is the measurement used by the controller.

## Incremental PID Control

The temperature error is:

$$
e(n) = T_{\mathrm{set}} - T_{\mathrm{avg}}(n)
$$

The controller uses the incremental PID form documented in the source:

$$
\Delta P(n) = K_p \left[e(n)-e(n-1)\right] + K_i e(n) + K_d \left[e(n)-2e(n-1)+e(n-2)\right]
$$

$$
P(n) = P(n-1) + \Delta P(n)
$$

This structure updates heater power directly from error increments instead of recomputing an absolute PID output each time.

## PID Gain Scheduling by Temperature Error

The code selects one of four gain zones using the absolute error \(|e|\), together with hysteresis:

- Zone 0: coarse control for \(|e| > 100\) K
- Zone 1: transition control for \(20 < |e| \le 100\) K
- Zone 2: fine control for \(5 < |e| \le 20\) K
- Zone 3: precision control for \(|e| \le 5\) K

Each zone has its own \((K_p, K_i, K_d)\) tuple and update frequency. This gives stronger action far from the setpoint and finer action near convergence.

## D-Term Low-Pass Filtering

The raw derivative contribution is filtered as:

$$
D_f(n) = \alpha_D D_{\mathrm{raw}}(n) + (1-\alpha_D) D_f(n-1)
$$

where \(\alpha_D = 0.3\) in the current code. This reduces oscillatory or noisy derivative response.

## Heater Power Split

The controller computes a total heater power \(P\), then splits it with a fixed ratio:

$$
P_{\mathrm{main}} = \mathrm{RATIO\_MAIN} \cdot P
$$

$$
P_{\mathrm{bottom}} = \left(1-\mathrm{RATIO\_MAIN}\right) \cdot P
$$

In the current implementation, `RATIO_MAIN = 0.8`, so the main heater receives 80% of the total power and the bottom heater receives 20%.

## UDM-Based Heat-Source Writing

The heater source terms do not compute power internally. Instead, `heater_pid_control` writes volumetric source density into `UDM[0]`:

$$
\dot{q}_{\mathrm{main}} = \frac{P_{\mathrm{main}}}{V_{\mathrm{main}}}
$$

$$
\dot{q}_{\mathrm{bottom}} = \frac{P_{\mathrm{bottom}}}{V_{\mathrm{bottom}}}
$$

The energy source hooks simply return `C_UDMI(c,t,UDM_QDOT)`, which keeps the control update and the source application decoupled.

## Crystal Rotation Momentum Source

When rotation is enabled, the UDF applies momentum sources in the silicon zone using the target rigid-body rotation:

$$
v_{\mathrm{rot},x} = -\omega z
$$

$$
v_{\mathrm{rot},z} = \omega x
$$

The source form is:

$$
S = -C_{\mathrm{rot}} \frac{(1-f_l)^2}{f_l^3 + \varepsilon} \left(v - v_{\mathrm{rot}}\right)
$$

where:

- \(f_l\) is the liquid fraction from `C_VOF`
- \(C_{\mathrm{rot}}\) is the Carman-Kozeny-style penalty coefficient
- \(\varepsilon\) prevents division by zero

The UDF applies this in both X and Z momentum equations.

## Carman-Kozeny Type Source Term

The factor

$$
\frac{(1-f_l)^2}{f_l^3 + \varepsilon}
$$

behaves like a porous-resistance penalty. It becomes small in liquid regions and large in solid-like regions, which drives the solid crystal toward the prescribed rotation while allowing the liquid to remain comparatively unconstrained.

## Liquid, Solid, and Mushy-Zone Behavior

- \(f_l = 1\): source tends to zero, so the liquid is not forced to rotate rigidly.
- \(0 < f_l < 1\): intermediate resistance creates a smooth mushy-zone transition.
- \(f_l = 0\): the penalty is large, so the solid tends toward the target rotational velocity.

This is why the source can be used in a single mixed silicon zone containing melt, mushy, and solid regions.

## Rotation Enable/Suspend State Machine

Rotation is not always active. The code uses a three-state controller:

- `OFF (0)` -> `ON (1)` after \(|e| < 10\) K for 50 consecutive PID updates
- `ON (1)` -> `SUSPENDED (2)` when \(|e| > 50\) K
- `SUSPENDED (2)` -> `ON (1)` after re-satisfying the same stability condition

This delays rotation until thermal control is sufficiently stable and suspends it if the interface temperature departs too far from the target.

## Fluent Parallel Reduction

Interface temperature sampling is parallel-compatible. Each process accumulates local sums, then global reductions are performed with `PRF_GRSUM1` for:

- temperature-area sum
- total selected area
- selected face count

That produces a consistent global area-weighted temperature without explicit host-to-node broadcasting.

## Persistent State File

The runtime state is written to `pid_state.pid`. In version `v3.0`, the file stores:

- total heater power
- previous PID errors
- filtered derivative history
- PID update count
- previous gain zone
- short temperature history
- current absolute error
- rotation state
- rotation stability counter

On restart, `heater_pid_control` tries to restore this state on its first call. This supports simulation continuation without resetting the controller.
