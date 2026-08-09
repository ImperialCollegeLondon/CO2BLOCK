use num_traits::Zero;
use std::num::NonZeroU64;
use uom::si::{
    angle::degree,
    compressibility_coefficient::per_megapascal,
    f64::*,
    length::kilometer,
    pressure::megapascal,
    ratio::{part_per_million, ratio},
    temperature_gradient::kelvin_per_kilometer,
};

use crate::calc::PressureGradient;

pub struct Args {
    pub reservoir: ReservoirParams,
    pub injection: InjectionParams,
    pub well_placement: WellPlaceParams,
    pub correction: Correction,
}

// TODO: implement partial defaults
pub struct ReservoirParams {
    pub domain: DomainParams,
    pub rock: RockParams,
    pub water: WaterParams,
    pub gas: GasParams,
}

pub struct WaterParams {
    pub compress: CompressibilityCoefficient,
    pub density: MassDensity,
    pub visc: DynamicViscosity,
}

pub struct GasParams {
    pub density: MassDensity,
    pub visc: DynamicViscosity,
}

pub struct RockParams {
    pub porosity: Ratio,
    pub permeability: Area,
    pub compress_rock: CompressibilityCoefficient,
    pub cohesion: Pressure,
    pub stress_ratio: Ratio,
    pub friction_angle: Angle,
    pub tensile_strength: Pressure,
    pub max_principal_stress: Pressure,
    pub top_pressure: Pressure,
}

pub struct DomainParams {
    pub thickness: Length,
    pub area: Area,
    pub domain_type: DomainType,
}

pub enum DomainType {
    Open,
    Closed,
}

pub enum Correction {
    Off,
    On,
}

pub struct InjectionParams {
    pub well_radius: Length,
    pub duration: Time,
    pub max_well_rate: MassRate,
}

pub struct WellPlaceParams {
    pub num_distance_samples: u64,
    pub inter_well_dist_min: Length,
    pub inter_well_dist_max: Option<Length>,
    pub num_wells_max: Option<NonZeroU64>,
}

pub trait DefaultProps {
    fn litho_grad() -> PressureGradient {
        Pressure::new::<megapascal>(23.) / Length::new::<kilometer>(1.)
    }

    fn hydro_grad() -> PressureGradient {
        Pressure::new::<megapascal>(10.) / Length::new::<kilometer>(1.)
    }

    fn temperature_grad() -> TemperatureGradient {
        TemperatureGradient::new::<kelvin_per_kilometer>(33.)
    }

    fn stress_ratio() -> Ratio {
        Ratio::new::<ratio>(0.7)
    }

    fn rock_friction_angle() -> Angle {
        Angle::new::<degree>(30.)
    }

    fn rock_cohesion() -> Pressure {
        Zero::zero()
    }

    fn rock_compressibility() -> CompressibilityCoefficient {
        CompressibilityCoefficient::new::<per_megapascal>(5e-4)
    }

    fn water_compressibility() -> CompressibilityCoefficient {
        CompressibilityCoefficient::new::<per_megapascal>(3e-4)
    }

    fn salinity() -> Ratio {
        Ratio::new::<part_per_million>(18e4)
    }
}
