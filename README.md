# CZ Silicon Furnace Fluent UDF

This repository contains a single ANSYS Fluent User-Defined Function (UDF) for a Czochralski (CZ) single-crystal silicon furnace CFD simulation. The UDF combines PID-based heater power control with a crystal-rotation momentum source and persistent runtime state handling.

## Main Features

- Incremental PID control of total heater power using interface temperature feedback
- Area-weighted interface temperature sampling over a specified radial band
- Gain scheduling based on temperature-error magnitude
- Low-pass-filtered derivative term for smoother PID response
- Fixed split of total heater power between main and bottom heaters
- UDM-based heat-source storage for Fluent energy source hooks
- Crystal rotation imposed through Carman-Kozeny-type X/Z momentum sources
- Rotation enable/suspend logic tied to PID convergence behavior
- Persistent PID and rotation state via `pid_state.pid`
- Parallel-compatible temperature reduction and Fluent console diagnostics

## Physical Background

The UDF targets a CZ silicon furnace model in which interface temperature is regulated toward a prescribed setpoint. Total heater power is adjusted by PID control, then distributed between the main heater and the bottom heater. A separate momentum source can impose crystal rotation in the silicon zone, with the source strength modulated by liquid fraction so the solid region follows the target rotation while the liquid region remains essentially free of rotational forcing.

## Core Workflow

1. Fluent evaluates `heater_pid_control` at the end of the iteration sequence.
2. The UDF samples interface-face temperature over the configured radial interval.
3. Parallel reductions form the global area-weighted average temperature.
4. The incremental PID law updates total heater power at the configured control frequency.
5. The power is split between heater zones and written to `UDM[0]`.
6. Energy source hooks read `UDM[0]` for the main and bottom heater zones.
7. Rotation state is updated from the temperature-error history.
8. Momentum source hooks apply crystal rotation only when rotation is enabled.
9. PID and rotation state are saved to `pid_state.pid` for resume behavior.

## Repository Structure

```text
cz-silicon-furnace-udf/
|-- README.md
|-- LICENSE
|-- .gitignore
|-- src/
|   `-- cz_pid_rotation_udf.c
|-- docs/
|   |-- setup.md
|   |-- theory.md
|   `-- changelog.md
`-- examples/
    `-- hook_configuration.md
```

## Requirements

- ANSYS Fluent with UDF compilation support
- C compilation environment compatible with Fluent on your platform
- One User-Defined Memory slot enabled in Fluent
- A case setup whose zone IDs and geometry assumptions match the UDF

The source header states it was written for Fluent 24.2, double precision (`3ddp`), Windows, and parallel execution. This repository does not claim compilation, runtime verification, or physical validation beyond the source code itself.

## Fluent Setup Summary

- Enable `1` User-Defined Memory slot before loading the UDF.
- Compile and load `src/cz_pid_rotation_udf.c`.
- Hook `heater_pid_control` as `Execute at End`.
- Hook `heater_main_source` to the main-heater energy source.
- Hook `heater_bottom_source` to the bottom-heater energy source.
- Hook `crystal_xmom_source` and `crystal_zmom_source` to the silicon-zone X/Z momentum sources.
- Hook `reset_pid` as `Execute On Demand` if you want a manual reset action.

See [`docs/setup.md`](docs/setup.md) and [`examples/hook_configuration.md`](examples/hook_configuration.md) for the detailed procedure.

## Important Parameters

Review these code parameters before use because they are case-specific:

| Parameter | Meaning |
| --- | --- |
| `ID_INTERFACE` | Interface face zone used for temperature feedback |
| `ID_HEATER_MAIN` | Main heater cell zone ID |
| `ID_HEATER_BOTTOM` | Bottom heater cell zone ID |
| `ID_SILICON` | Silicon melt + crystal cell zone ID |
| `V_HEATER_MAIN` | Main heater zone volume used for volumetric source density |
| `V_HEATER_BOTTOM` | Bottom heater zone volume used for volumetric source density |
| `R_MIN` / `R_MAX` | Radial window for interface averaging |
| `T_SET` | Interface temperature target |
| `ROT_OMEGA` | Target crystal angular speed |

## Notes and Limitations

- The repository preserves the original UDF logic and comments; only project organization and documentation were added.
- Zone IDs, physical parameters, PID equations, and rotation-source equations were not altered.
- The UDF assumes one silicon zone for the melt/crystal region when applying the rotation source.
- The UDF writes and reads `pid_state.pid` in the Fluent working directory during runtime.
- No Fluent case, data, mesh, figure, benchmark, or validation artifact is included here.

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).
