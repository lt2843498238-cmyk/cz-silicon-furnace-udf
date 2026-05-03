# Setup Guide

## Environment Assumptions

- The UDF is written for ANSYS Fluent with compiled UDF support.
- The source header identifies the original target environment as Fluent 24.2, double precision (`3ddp`), Windows, and parallel execution.
- The simulation is expected to contain:
  - an interface face zone for temperature feedback
  - a main heater cell zone
  - a bottom heater cell zone
  - a silicon melt + crystal cell zone for rotation sources

## Parameters to Check in the Code

Before using the UDF, verify these definitions in [`src/cz_pid_rotation_udf.c`](../src/cz_pid_rotation_udf.c):

```c
ID_INTERFACE
ID_HEATER_MAIN
ID_HEATER_BOTTOM
ID_SILICON
V_HEATER_MAIN
V_HEATER_BOTTOM
R_MIN
R_MAX
T_SET
ROT_OMEGA
```

These values must match your mesh topology, heater volumes, and intended operating point.

## Required User-Defined Memory Setting

The UDF requires exactly one User-Defined Memory slot:

1. In Fluent, open `Define -> User-Defined -> Memory`.
2. Set the number of UDMs to `1`.
3. Apply the setting before compiling/loading the UDF.

`UDM[0]` is used to store volumetric heat-source density for the heater source functions.

## Compile and Load the UDF

1. Open Fluent and load your case.
2. Enable one UDM slot as described above.
3. Open the compiled UDF panel.
4. Add `src/cz_pid_rotation_udf.c` to the source list.
5. Build the library with your Fluent-compatible compiler.
6. Load the compiled library into the case.

This repository does not claim that compilation has been performed successfully in your environment.

## Function Hooks

Hook the UDF functions as follows:

- `Execute at End` -> `heater_pid_control`
- `Energy source for main heater zone` -> `heater_main_source`
- `Energy source for bottom heater zone` -> `heater_bottom_source`
- `X-Momentum source for silicon zone` -> `crystal_xmom_source`
- `Z-Momentum source for silicon zone` -> `crystal_zmom_source`
- `Execute On Demand` -> `reset_pid`

An explicit hook table is also provided in [`examples/hook_configuration.md`](../examples/hook_configuration.md).

## Zones That Need Source Terms

- Main heater zone: energy source `heater_main_source`
- Bottom heater zone: energy source `heater_bottom_source`
- Silicon melt + crystal zone: X-momentum `crystal_xmom_source`
- Silicon melt + crystal zone: Z-momentum `crystal_zmom_source`

The source code currently defines these zone IDs:

- `ID_INTERFACE = 4`
- `ID_HEATER_MAIN = 4508`
- `ID_HEATER_BOTTOM = 4505`
- `ID_SILICON = 4592`

## How to Resume a Simulation

The UDF stores PID and rotation state in `pid_state.pid`.

To resume:

1. Open the existing case and data files.
2. Load the compiled UDF library again.
3. Keep `pid_state.pid` in the Fluent working directory.
4. Continue the calculation.

On the first `Execute at End` call, the UDF attempts to restore the saved state and re-write the current heater source values to `UDM[0]`.

## How to Reset PID State

Use the on-demand function `reset_pid` when you want to clear the persisted control history without editing the C source.

Effects of `reset_pid`:

- resets total power to `P_START`
- clears PID error history and filtered derivative state
- clears temperature history
- turns rotation off
- resets the rotation stability counter
- writes fresh heater source values into the heater-zone UDM

For a fresh start from disk, remove `pid_state.pid` before loading/running the case.

## Common Mistakes and Troubleshooting

- `No UDM defined`: the heater source functions expect `UDM[0]`; define one UDM before loading the UDF.
- `No interface averaging`: if the code prints a warning about zero face area, check `ID_INTERFACE`, `R_MIN`, and `R_MAX`.
- `Wrong heaters receive power`: confirm `ID_HEATER_MAIN`, `ID_HEATER_BOTTOM`, `V_HEATER_MAIN`, and `V_HEATER_BOTTOM`.
- `Rotation source has no effect`: confirm the momentum hooks are assigned to the silicon zone and that the rotation state has actually turned on.
- `Rotation never enables`: verify that the interface error can remain below the enable threshold long enough for the consecutive-count logic.
- `Unexpected restart behavior`: if old controller state is being restored, remove `pid_state.pid` or call `reset_pid`.
- `Compilation issues`: check that your Fluent version, precision mode, and compiler toolchain are consistent with compiled UDF requirements.
