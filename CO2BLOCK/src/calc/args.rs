use std::num::NonZeroU64;
use uom::si::f64::*;

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
