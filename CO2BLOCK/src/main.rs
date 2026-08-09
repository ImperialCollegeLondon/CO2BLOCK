use co2block::calc::InputReservoirParams;
use uom::si::{
    f64::{MassDensity, Pressure},
    mass_density::kilogram_per_cubic_meter,
    pressure::bar,
};

fn main() {
    let mut data = InputReservoirParams::default();
    data.rock_compress = Some(co2block::calc::default_props::rock_compressibility());
    data.pressure_top = Some(Pressure::new::<bar>(10.));
    data.co2_density = Some(MassDensity::new::<kilogram_per_cubic_meter>(1.));
    let json = serde_json::to_string_pretty(&data).unwrap();
    println!("{json}");
    {
        let mut wtr = csv::Writer::from_path("input-template.csv").unwrap();
        wtr.serialize(data).unwrap();
    }
    let mut rdr = csv::Reader::from_path("input-template.csv").unwrap();
    let _data: InputReservoirParams = rdr.deserialize().into_iter().next().unwrap().unwrap();
}
