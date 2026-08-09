use std::{io, str::FromStr};

use co2block::calc::InputReservoirParams;
use uom::si::f64::Length;

fn main() {
    let data = InputReservoirParams::default();
    let json = serde_json::to_string_pretty(&data).unwrap();
    println!("{json}");
    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(data).unwrap()
}
