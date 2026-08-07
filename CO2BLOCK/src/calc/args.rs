use std::num::NonZeroU64;
use uom::si::f64::*;

pub struct _Args {
    pub reservoir: ReservoirParams,
    pub injection: InjectionParams,
    well_placement: _WellPlaceParams,
    pub correction: Correction,
}

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
    pub duration_injection: Time,
    pub rate_max: MassRate,
}

struct _WellPlaceParams {
    pub num_distance_samples: u64,
    pub inter_well_dist_min: Length,
    pub inter_well_dist_max: Option<Length>,
    pub num_wells_max: Option<NonZeroU64>,
}

uom::unit! {
    system: uom::si;
    quantity: uom::si::mass_rate;
    @megaton_per_year: prefix!(mega) * 1.0_E3 / 3.1536_E7; "Mt/y", "megaton per year", "megatons per year";
}
