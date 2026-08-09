use std::{io, str::FromStr};

use co2block::calc::InputReservoirParams;

fn main() {
    let mut rdr = csv::Reader::from_path("input-template.csv").unwrap();
    let data: InputReservoirParams = rdr.deserialize().into_iter().next().unwrap().unwrap();
    let json = serde_json::to_string_pretty(&data).unwrap();
    println!("{json}");
    let mut wtr = csv::Writer::from_path("input-template.csv").unwrap();
    wtr.serialize(data).unwrap()
}
