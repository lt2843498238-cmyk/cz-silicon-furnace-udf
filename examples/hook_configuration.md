# Example Hook Configuration

Use the following Fluent hook mapping as the baseline configuration for this UDF:

| Fluent location                | UDF function           |
| ------------------------------ | ---------------------- |
| Execute at End                 | `heater_pid_control`   |
| Main heater Energy source      | `heater_main_source`   |
| Bottom heater Energy source    | `heater_bottom_source` |
| Silicon zone X-Momentum source | `crystal_xmom_source`  |
| Silicon zone Z-Momentum source | `crystal_zmom_source`  |
| Execute On Demand              | `reset_pid`            |

Notes:

- The main heater, bottom heater, and silicon zone IDs must match the values defined in the source file.
- The interface feedback zone is referenced internally by `ID_INTERFACE`; it is not a separate hook assignment.
