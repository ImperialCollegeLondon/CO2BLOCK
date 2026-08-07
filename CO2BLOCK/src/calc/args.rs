use std::num::NonZeroU64;
use uom::si::f64::*;

pub struct Args {
    pub hyper_params: HyperParams,
    pub injection_params: InjectionParams,
    pub domain_params: DomainParams,
    pub rock_params: RockParams,
    pub water_params: WaterParams,
    pub gas_params: GasParams,
    pub well_placement: WellPlaceParams,
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

pub struct HyperParams {
    pub correction: Correction,
    pub num_distance_samples: u64,
}

pub struct InjectionParams {
    pub well_radius: Length,
    pub duration_injection: Time,
    pub rate_max: MassRate,
}

pub struct WellPlaceParams {
    pub inter_well_dist_min: Length,
    pub inter_well_dist_max: Option<Length>,
    pub num_wells_max: Option<NonZeroU64>,
}

uom::unit! {
    system: uom::si;
    quantity: uom::si::mass_rate;
    @megaton_per_year: {const{prefix!(mega) * 1.0_E3 / 3.1536_E7}}; "Mt/y", "megaton per year", "megatons per year";
}
