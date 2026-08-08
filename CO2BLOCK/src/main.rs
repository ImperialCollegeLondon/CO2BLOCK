use co2block::calc::*;
use uom::si::{
    area::square_kilometer,
    compressibility_coefficient::per_megapascal,
    f64::*,
    length::{kilometer, meter},
    pressure::megapascal,
    ratio::ratio,
    time::year,
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
        injection: InjectionParams {
            well_radius: Length::new::<meter>(0.2),
            duration: Time::new::<year>(100.),
            max_well_rate: MassRate::new::<megaton_per_year>(20.),
        },
        well_placement: WellPlaceParams {
            num_distance_samples: 0,
            inter_well_dist_min: Length::new::<kilometer>(2.),
            inter_well_dist_max: None,
            num_wells_max: None,
        },
        correction: Correction::Off,
    };
}
