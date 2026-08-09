use co2block::calc::InputReservoirParams;

fn main() {
    let mut data = InputReservoirParams::default();
    data.rock_compress = Some(co2block::calc::default_props::rock_compressibility());
    let json = serde_json::to_string_pretty(&data).unwrap();
    println!("{json}");
    {
        let mut wtr = csv::Writer::from_path("input-template.csv").unwrap();
        wtr.serialize(data).unwrap();
    }
    let mut rdr = csv::Reader::from_path("input-template.csv").unwrap();
    let _data: InputReservoirParams = rdr.deserialize().into_iter().next().unwrap().unwrap();
}
