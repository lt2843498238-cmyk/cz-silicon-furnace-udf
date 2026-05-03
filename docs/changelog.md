# Changelog

## v3.0

- Added PID-based heater power control with incremental update logic.
- Added crystal rotation momentum sources for the silicon zone.
- Added persistent runtime state through `pid_state.pid`.
- Added parallel-compatible interface temperature reduction using Fluent reduction macros.
- Added structured diagnostic console output for temperature, power, PID terms, and rotation status.
