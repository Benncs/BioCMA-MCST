\page kernels List of Kokkkos Kernels


# List of Kokkos Kernels

## Particle Container Management (MC)
The following kernels manage the particle container, including initialization, repartitioning, merging, and defragmentation.

| Type         | Name                  | Policy                                      | Brief Description                                                                 | Number of Calls          |
|--------------|-----------------------|---------------------------------------------|-----------------------------------------------------------------------------------|--------------------------|
| reduction    | `mc_init_first`       | `range(size)`                               | Spawns particles by setting positions, calling model initialization, and calculating total mass. | 1                        |
| for          | `get_repartition`     | `range(size)`                               | Counts the number of particles per compartment (if multiple compartments exist).   | `n_export`               |
| for          | `insert_merge`        | `TeamPolicy(n_add_item, Kokkos::AUTO, Model::n_var)` | Back-inserts to merge the main container and buffer.                              | `n_step`                 |
| scan         | `find_and_fill_gap`   | `range(size)`                               | Finds non-idle particles and replaces them with new ones (defragmentation).      | If `n_non_idle > threshold` |

---

## Iteration Cycle Kernels
Each iteration cycle launches the following kernels:

| Type         | Name                  | Policy                                      | Brief Description                                                                 | Number of Calls          |
|--------------|-----------------------|---------------------------------------------|-----------------------------------------------------------------------------------|--------------------------|
| for          | `cycle_move`          | `team: (static_cast<int>(range) + team_size - 1)/team_size, team_size` | Moves particles based on the flowmap (if `n_compartment > 1`).                   | `n_step`                 |
| reduce       | `cycle_move_leave`    | `team: (256, auto)`                         | Returns the number of particles leaving (if continuous reactor with `feed != 0`). | `n_step`                 |
| reduce       | `cycle_model`         | `team: (static_cast<int>(range) + team_size - 1)/team_size, team_size` | Updates the model, handles division, and returns the number of particles leaving and waiting for allocation. | `n_step`                 |
| reduce       | `cycle_scatter`       | `team: range(size)`                         | Scatters particle contributions.                                                 | `n_step`                 |

---

## Data Export Kernels
The following kernel is used for exporting data:

| Type         | Name                  | Policy                                      | Brief Description                                                                 | Number of Calls          |
|--------------|-----------------------|---------------------------------------------|-----------------------------------------------------------------------------------|--------------------------|
| for          | `get_properties`      | `range(size)`                               | Deep copies the particle container to the host in a view format for HDF5 export. Skips non-idle particles. | n_export               |
