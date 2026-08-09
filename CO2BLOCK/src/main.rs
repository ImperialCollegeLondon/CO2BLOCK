use co2block::calc::{InputReservoirParams, ReservoirParams};

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
