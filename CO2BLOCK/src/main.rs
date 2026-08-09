use std::{io, str::FromStr};

use co2block::calc::InputReservoirParams;

fn main() {
    let data = InputReservoirParams::default();
    let json = serde_json::to_string_pretty(&data).unwrap();
    println!("{json}");
    let mut wtr = csv::Writer::from_path("input-template.csv").unwrap();
    wtr.serialize(data).unwrap()
}
