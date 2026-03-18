use std::env;
use std::io::{self, IsTerminal, Write};
use std::time::{Duration, Instant};

pub struct ProgressReporter {
    prefix: String,
    enabled: bool,
}

impl ProgressReporter {
    pub fn new(prefix: impl Into<String>) -> Self {
        let enabled = io::stderr().is_terminal()
            && env::var_os("KAROSPACE_COMPANION_NO_PROGRESS").is_none();
        Self {
            prefix: prefix.into(),
            enabled,
        }
    }

    pub fn stage<'a>(&'a self, name: impl Into<String>) -> ProgressStage<'a> {
        let name = name.into();
        if self.enabled {
            eprintln!("[{}] {}", self.prefix, name);
        }
        ProgressStage {
            reporter: self,
            name,
            started_at: Instant::now(),
            live_progress: false,
            finished: false,
        }
    }
}

pub struct ProgressStage<'a> {
    reporter: &'a ProgressReporter,
    name: String,
    started_at: Instant,
    live_progress: bool,
    finished: bool,
}

impl ProgressStage<'_> {
    pub fn note(&mut self, message: impl AsRef<str>) {
        if !self.reporter.enabled {
            return;
        }
        self.clear_live_progress();
        eprintln!(
            "[{}] {}: {}",
            self.reporter.prefix,
            self.name,
            message.as_ref()
        );
    }

    pub fn progress(&mut self, label: impl AsRef<str>, current: usize, total: usize) {
        if !self.reporter.enabled || total <= 1 || current == 0 {
            return;
        }
        if !should_emit_progress(current, total) {
            return;
        }
        eprint!(
            "\r[{}] {}: {} {}/{}",
            self.reporter.prefix,
            self.name,
            label.as_ref(),
            current,
            total
        );
        let _ = io::stderr().flush();
        self.live_progress = true;
    }

    pub fn finish(mut self, detail: impl AsRef<str>) {
        if self.reporter.enabled {
            self.clear_live_progress();
            let detail = detail.as_ref();
            let elapsed = format_elapsed(self.started_at.elapsed());
            if detail.is_empty() {
                eprintln!(
                    "[{}] {} done in {}",
                    self.reporter.prefix, self.name, elapsed
                );
            } else {
                eprintln!(
                    "[{}] {} done in {} ({})",
                    self.reporter.prefix, self.name, elapsed, detail
                );
            }
        }
        self.finished = true;
    }

    fn clear_live_progress(&mut self) {
        if self.live_progress {
            eprintln!();
            self.live_progress = false;
        }
    }
}

impl Drop for ProgressStage<'_> {
    fn drop(&mut self) {
        if !self.reporter.enabled || self.finished {
            return;
        }
        self.clear_live_progress();
        eprintln!(
            "[{}] {} stopped after {}",
            self.reporter.prefix,
            self.name,
            format_elapsed(self.started_at.elapsed())
        );
    }
}

fn should_emit_progress(current: usize, total: usize) -> bool {
    if current == total || current == 1 {
        return true;
    }
    let step = (total / 20).max(1);
    current % step == 0
}

fn format_elapsed(duration: Duration) -> String {
    let secs = duration.as_secs();
    if secs >= 60 {
        format!("{}m{:02}s", secs / 60, secs % 60)
    } else if secs > 0 {
        format!("{}.{:01}s", secs, duration.subsec_millis() / 100)
    } else {
        format!("{}ms", duration.as_millis())
    }
}
