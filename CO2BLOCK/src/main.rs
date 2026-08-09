use std::str::FromStr;

use co2block::calc::{InjectionParams, InputReservoirParams, ReservoirParams};
use uom::si::f64::{Length, MassRate, Time};

fn main() -> anyhow::Result<()> {
    let input: InputReservoirParams = csv::Reader::from_path("input-template.csv")?
        .deserialize()
        .into_iter()
        .next()
        .ok_or(anyhow::anyhow!("Empty input"))??;
    let reservoir: ReservoirParams = input.into();
    println!("{reservoir:#?}");

    let step = Length::from_str("1 km")?;

    let injection = InjectionParams {
        well_radius: Length::from_str("0.2 m")?,
        duration: Time::from_str("100 y")?,
        max_well_rate: MassRate::from_str("20 Mt/y")?,
    };

    co2block::calc::co2block_with_placement(
        step,
        1,
        reservoir,
        injection,
        co2block::calc::Correction::Off,
    );

    Ok(())
}
