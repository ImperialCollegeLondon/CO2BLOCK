use std::str::FromStr;

use co2block::calc::{InjectionParams, InputReservoirParams, ReservoirParams, megaton_per_year};
use uom::si::{
    f64::{Length, MassRate, Time},
    length::{kilometer, meter},
    time::year,
};

fn main() -> anyhow::Result<()> {
    let input: InputReservoirParams = csv::Reader::from_path("input-template.csv")?
        .deserialize()
        .into_iter()
        .next()
        .ok_or(anyhow::anyhow!("Empty input"))??;
    let reservoir: ReservoirParams = input.into();
    // println!("{reservoir:#?}");

    let step = Length::new::<kilometer>(2.);

    let injection = InjectionParams {
        well_radius: Length::new::<meter>(0.2),
        duration: Time::new::<year>(100.),
        max_well_rate: MassRate::new::<megaton_per_year>(20.),
    };

    let well_rate = co2block::calc::co2block_with_placement(
        step,
        1,
        reservoir,
        injection,
        co2block::calc::Correction::Off,
    );
    println!(
        "{}",
        well_rate
            .into_format_args(megaton_per_year, uom::fmt::DisplayStyle::Abbreviation)
            .to_string()
    );

    Ok(())
}
