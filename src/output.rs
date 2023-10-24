use std::io::{stdout, Error, Write};

use crate::{cli::Config, process::ProcessResults};

pub fn output_results(cfg: &Config, res: &ProcessResults) -> Result<(), Error> {
    let mut wrt = stdout();
    writeln!(&mut wrt, "{{")?;
    let in_list = cfg.fli().json_output(&mut wrt, 3, false)?;
    res.json_output(&mut wrt, cfg.trim(), 3, in_list)?;
    writeln!(&mut wrt, "\n}}")
}
