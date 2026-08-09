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
    thermodynamic_temperature::degree_celsius,
};

use crate::calc::{
    PressureGradient,
    eos::{self, gas_density, gas_viscosity},
};

// TODO implement serde and clap for hyper-parameters
pub struct Args {
    pub reservoir: ReservoirParams,
    pub injection: InjectionParams,
    pub well_placement: WellPlaceParams,
    pub correction: Correction,
    pub boundary: Boundary,
}

pub enum Boundary {
    Closed,
    Open,
}

pub struct ReservoirParams {
    pub domain: DomainParams,
    pub rock: RockParams,
    pub water: WaterParams,
    pub gas: GasParams,
}

pub struct WaterParams {
    pub compress: CompressibilityCoefficient,
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

mod default_props {
    use super::*;
    pub(crate) fn litho_grad() -> PressureGradient {
        Pressure::new::<megapascal>(23.) / Length::new::<kilometer>(1.)
    }

    pub(crate) fn hydro_grad() -> PressureGradient {
        Pressure::new::<megapascal>(10.) / Length::new::<kilometer>(1.)
    }

    pub(crate) fn temperature_grad() -> TemperatureGradient {
        TemperatureGradient::new::<kelvin_per_kilometer>(33.)
    }

    pub(crate) fn stress_ratio() -> Ratio {
        Ratio::new::<ratio>(0.7)
    }

    pub(crate) fn rock_friction_angle() -> Angle {
        Angle::new::<degree>(30.)
    }

    pub(crate) fn rock_cohesion() -> Pressure {
        Zero::zero()
    }

    pub(crate) fn rock_compressibility() -> CompressibilityCoefficient {
        CompressibilityCoefficient::new::<per_megapascal>(5e-4)
    }

    pub(crate) fn water_compressibility() -> CompressibilityCoefficient {
        CompressibilityCoefficient::new::<per_megapascal>(3e-4)
    }

    pub(crate) fn salinity() -> Ratio {
        Ratio::new::<part_per_million>(18e4)
    }
}

// TODO: implement serde from formatted strings
// TODO: implement serde-csv
pub struct InputReservoirParams {
    pub shallowest_depth: Length,
    pub mean_depth: Length,
    pub thickness: Length,
    pub area: Area,
    pub permeability: Area,
    pub porosity: Ratio,
    pub rock_compress: Option<CompressibilityCoefficient>,
    pub water_compress: Option<CompressibilityCoefficient>,
    pub pressure_top: Option<Pressure>,
    pub pressure_center: Option<Pressure>,
    pub co2_density: Option<MassDensity>,
    pub co2_viscosity: Option<DynamicViscosity>,
    pub water_viscosity: Option<DynamicViscosity>,
    pub salinity: Option<Ratio>,
    pub temperature_center: Option<ThermodynamicTemperature>,
    pub total_max_princ_stress: Option<Pressure>,
    pub princ_stress_ratio: Option<Ratio>,
    pub rock_friction_angle: Option<Angle>,
    pub rock_cohesion: Option<Pressure>,
    pub rock_tensile_strength: Option<Pressure>,
}

impl From<InputReservoirParams> for ReservoirParams {
    fn from(input: InputReservoirParams) -> Self {
        let domain = DomainParams {
            thickness: input.thickness,
            area: input.area,
        };
        let cohesion = input
            .rock_cohesion
            .unwrap_or(default_props::rock_cohesion());
        let rock = RockParams {
            porosity: input.porosity,
            permeability: input.permeability,
            compress_rock: input
                .rock_compress
                .unwrap_or(default_props::rock_compressibility()),
            cohesion,
            stress_ratio: input
                .princ_stress_ratio
                .unwrap_or(default_props::stress_ratio()),
            friction_angle: input
                .rock_friction_angle
                .unwrap_or(default_props::rock_friction_angle()),
            tensile_strength: input.rock_tensile_strength.unwrap_or(cohesion / 2.),
            max_principal_stress: input
                .total_max_princ_stress
                .unwrap_or(default_props::litho_grad() * input.mean_depth),
            top_pressure: input
                .pressure_top
                .unwrap_or(default_props::hydro_grad() * input.mean_depth),
        };

        let temperature = input.temperature_center.unwrap_or(
            ThermodynamicTemperature::new::<degree_celsius>(15.)
                + default_props::temperature_grad() * input.mean_depth,
        );

        let salinity = input.salinity.unwrap_or(default_props::salinity());

        let water = WaterParams {
            compress: input
                .water_compress
                .unwrap_or(default_props::water_compressibility()),
            visc: input
                .water_viscosity
                .unwrap_or(eos::brine_viscosity(temperature, salinity)),
        };

        let pressure = input
            .pressure_center
            .unwrap_or(default_props::hydro_grad() * input.mean_depth);

        let gas_density = gas_density(pressure, temperature);
        let gas = GasParams {
            density: gas_density,
            visc: gas_viscosity(temperature, gas_density),
        };

        ReservoirParams {
            domain,
            rock,
            water,
            gas,
        }
    }
}
