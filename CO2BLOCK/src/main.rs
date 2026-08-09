use co2block::calc::{InputReservoirParams, ReservoirParams, default_props};
use uom::si::{
    dynamic_viscosity::{DynamicViscosity, pascal_second},
    f64::{Length, MassDensity, Pressure, ThermodynamicTemperature},
    length::meter,
    mass_density::kilogram_per_cubic_meter,
    pressure::{bar, megapascal},
    thermodynamic_temperature::degree_celsius,
};

fn main() -> anyhow::Result<()> {
    let input: InputReservoirParams = csv::Reader::from_path("input-template.csv")?
        .deserialize()
        .into_iter()
        .next()
        .ok_or(anyhow::anyhow!("Empty input"))??;
    let input_res: ReservoirParams = input.into();
    println!("{input_res:#?}");
    Ok(())
}
