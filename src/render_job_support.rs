pub fn resolve_render_max_jobs(cli: Option<usize>, cfg: Option<usize>, threads: usize) -> usize {
    cli.or(cfg).unwrap_or(threads.max(1)).max(1)
}
