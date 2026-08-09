use co2block::calc::InputReservoirParams;

use co2block::calc::milli_darcy;

fn main() {
    let data = InputReservoirParams::default();
    let json = serde_json::to_string_pretty(&data).unwrap();
    println!("{json}");
    {
        let mut wtr = csv::Writer::from_path("input-template.csv").unwrap();
        wtr.serialize(data).unwrap();
    }
    let mut rdr = csv::Reader::from_path("input-template.csv").unwrap();
    let _data: InputReservoirParams = rdr.deserialize().into_iter().next().unwrap().unwrap();
}
