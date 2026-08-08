use co2block::calc::*;
use uom::si::{
    area::square_kilometer, compressibility_coefficient::per_megapascal, f64::*, length::meter,
    pressure::megapascal, ratio::ratio,
};

fn main() {
    let args = co2block::calc::Args {
        reservoir: ReservoirParams {
            domain: DomainParams {
                thickness: Length::new::<meter>(100.),
                area: Area::new::<square_kilometer>(1_000.),
                domain_type: DomainType::Open,
            },
            rock: RockParams {
                porosity: Ratio::new::<ratio>(0.2),
                permeability: Area::new::<milli_darcy>(100.),
                compress_rock: CompressibilityCoefficient::new::<per_megapascal>(5e-4),
                max_principal_stress: Pressure::new::<megapascal>(50.),
                ..Default::default()
            },
            water: WaterParams {
                compress: CompressibilityCoefficient::new::<per_megapascal>(3e-4),
                ..Default::default()
            },
            gas: Default::default(),
        },
        injection: (),
        well_placement: (),
        correction: (),
    };
}
