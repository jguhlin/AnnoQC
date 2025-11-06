# Bevy Integration Notes

- We use Bevy 0.17’s minimal `App` (no default plugins) to orchestrate the analysis pipeline. `Startup` ingests jobs, `Update` dispatches DIAMOND/MAFFT work onto `AsyncComputeTaskPool`, and completion is signalled via a `CompletionState` resource.
- When adding new systems, keep them deterministic; avoid side effects outside resources because the CLI runs `app.update()` in a tight loop until `CompletionState.finished` flips to `true`.
- Long-running tasks should be spawned through `AsyncComputeTaskPool::spawn`. Remember to poll with `futures_lite::future::poll_once` and drop completed `Task`s (`swap_remove` + `let _ = ...`) so they do not cancel prematurely.
- If you bump Bevy, revisit `PipelineConfigResource` and the system registrations (`add_systems(Startup, …)` / `add_systems(Update, …)`)—the scheduling APIs changed rapidly pre-0.17.
- For in-flight limits, we rely on `max_in_flight = cfg.threads`. Increase the number of threads through the CLI/config instead of spawning extra systems.

## Common Pitfall (0.17.x)

- AsyncComputeTaskPool must be initialized before use. If you see:
  "The AsyncComputeTaskPool has not been initialized yet. Please call AsyncComputeTaskPool::get_or_init beforehand.",
  initialize it explicitly before building the `App`:

  `AsyncComputeTaskPool::get_or_init(|| TaskPoolBuilder::new().num_threads(cfg.threads).build());`

  Then spawn tasks via `AsyncComputeTaskPool::get().spawn(...)` in systems.
