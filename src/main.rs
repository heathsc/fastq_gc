#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

use anyhow::Context;
use crossbeam_channel::unbounded;
use crossbeam_utils::thread;

mod cli;
mod control_seq;
mod output;
mod process;
mod reader;
mod utils;

use process::{process_thread, ProcessResults};

fn main() -> anyhow::Result<()> {
    let cfg = cli::handle_cli()?;

    let nt = cfg.threads();

    let mut error = false;
    let mut res = ProcessResults::new(cfg.control_seq().map(|c| c.seq_ids()));

    thread::scope(|scope| {
        // Channel used to send full buffers to process threads
        let (full_send, full_recv) = unbounded();

        // Channel used to send and receive empty buffers
        let (empty_send, empty_recv) = unbounded();

        let mut process_tasks = Vec::with_capacity(nt);
        for ix in 0..nt {
            let rx = full_recv.clone();
            let tx = empty_send.clone();
            let cfg = &cfg;
            process_tasks.push(scope.spawn(move |_| process_thread(cfg, ix, rx, tx)));
        }

        drop(full_recv);
        drop(empty_send);

        if let Err(e) = reader::reader(&cfg, empty_recv, full_send) {
            error!("{:?}", e);
            error = true;
        }

        // Wait for analysis threads
        for jh in process_tasks.drain(..) {
            match jh.join().expect("Error joining analysis thread thread") {
                Err(e) => {
                    error!("{:?}", e);
                    error = true
                }
                Ok(r) => res += r,
            }
        }
    })
    .expect("Error in scope generation");

    if error {
        Err(anyhow!("Error occurred during processing"))
    } else {
        output::output_results(&cfg, &res)
            .with_context(|| "Error occurred while writing output JSON file")
    }
}
