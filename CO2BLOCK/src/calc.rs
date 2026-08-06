#![allow(unused)]

use std::num::NonZeroU64;

use uom::si::{
    length::meter,
    mass_rate::ton_per_day,
    ratio::ratio,
    time::{day, year},
};

#[allow(non_snake_case)]
pub fn nordbotten_impl(r: Length, R: Length, psi: Length, R_ext: Length, gamma: Ratio) -> Ratio {
    if (r > psi) && (r > R) {
        return Ratio::new::<uom::si::ratio::ratio>(0.);
    }

    if r <= R {
        return fd_nor(r, R, R_ext);
    }

    fd_nor(psi, R, R_ext) + gamma * (psi / r).ln()
}

#[allow(non_snake_case)]
fn fd_nor(x: Length, R: Length, R_ext: Length) -> Ratio {
    if x >= R_ext {
        return Ratio::new::<uom::si::ratio::ratio>(0.);
    }

    if R <= R_ext {
        return (R / x).ln();
    }

    (R_ext / x).ln()
        + Ratio::new::<uom::si::ratio::ratio>(2. / 2.25) * (R / R_ext).powi(uom::typenum::P2::new())
        - Ratio::new::<uom::si::ratio::ratio>(3. / 4.)
}

mod cached;

use uom::si::f64::*;

pub struct Args {
    pub correction: Correction,
    pub inter_well_dist_min: Len<kilometer>,
    pub inter_well_dist_max: Option<Len<kilometer>>,
    pub num_distances: NonZeroU64,
    pub num_wells_max: Option<NonZeroU64>,
    pub well_radius: Len<meter>,
    pub duration_injection: Time<year>,
    pub rate_max: MegatonsPerYear,
    pub thickness: Len<meter>,
    pub area: Area,         // 10^-15 m2
    pub permeability: Area, //km2
    pub porosity: Ratio,
    pub compress_rock: CompressibilityCoefficient, // 1/Pa (MPa in tables)
    pub compress_water: CompressibilityCoefficient, // 1/Pa
    pub stress_ratio: Ratio,
    pub rock_friction_angle: Angle,      // deg
    pub rock_tensile_strength: Pressure, // MPa
    pub rock_cohesion: Pressure,         // MPa
}

fn gamma(v_c: DynamicViscosity, v_w: DynamicViscosity) -> Ratio {
    v_c / v_w
}

fn delta(v_c: DynamicViscosity, v_w: DynamicViscosity) -> Ratio {
    Ratio::new::<ratio>(1.) - gamma(v_c, v_w)
}

pub enum DomainType {
    Open,
    Closed,
}

pub struct MegatonsPerYear {
    _value: uom::si::f64::MassRate,
}

uom::unit! {
    system: uom::si;
    quantity: uom::si::mass_rate;
    @megaton_per_year: {const{prefix!(mega) * 1.0_E3 / 3.1536_E7}}; "Mt/y", "megaton per year", "megatons per year";
}

impl MegatonsPerYear {
    pub fn new(value: f64) -> Self {
        MegatonsPerYear {
            _value: uom::si::f64::MassRate::new::<megaton_per_year>(value),
        }
    }
}

impl std::fmt::Display for MegatonsPerYear {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self._value
            .into_format_args(ton_per_day, uom::fmt::DisplayStyle::Abbreviation)
            .fmt(f)
    }
}

pub struct Time<Unit> {
    _value: uom::si::f64::Time,
    _unit: std::marker::PhantomData<Unit>,
}

impl<Unit: uom::si::time::Unit + uom::Conversion<f64, T = f64>> Time<Unit> {
    pub fn new(value: f64) -> Self {
        Time {
            _value: uom::si::f64::Time::new::<Unit>(value),
            _unit: std::marker::PhantomData,
        }
    }
}

pub struct Len<Unit> {
    _value: Length,
    _unit: std::marker::PhantomData<Unit>,
}

impl<Unit: uom::si::length::Unit + uom::Conversion<f64, T = f64>> Len<Unit> {
    pub fn new(value: f64) -> Self {
        Len {
            _value: Length::new::<Unit>(value),
            _unit: std::marker::PhantomData,
        }
    }
}

use uom::si::length::kilometer;

pub enum Correction {
    Off,
    On,
}

fn calc_total_compress(
    c_rock: CompressibilityCoefficient,
    poro: Ratio,
    c_water: CompressibilityCoefficient,
) -> CompressibilityCoefficient {
    c_rock + poro * c_water
}

fn calc_influence_radius(
    perm: Area,
    time: uom::si::f64::Time,
    visc_w: DynamicViscosity,
    compr_total: CompressibilityCoefficient,
) -> Length {
    (perm * time / (visc_w * compr_total) * 2.246).sqrt()
}

fn calc_max_num_wells(area: Area, dist_min: Length) -> uom::si::u32::Ratio {
    let res_float = (area / (dist_min * dist_min));
    use num_traits::ToPrimitive;
    let res_raw: u32 = res_float
        .value
        .to_u32()
        .expect("divison + floor result must be valid & round");
    uom::si::u32::Ratio::new::<ratio>(res_raw)
}

fn calc_well_dist_max(area: Area) -> Length {
    (area * 2.).sqrt() / 2.
}

pub fn calculate(args: Args) {
    // arrange wells in a ~square grid, making up to
    // sqrt(max_num) columns in total
    // num_x * (num_x +? 1)
}

mod tests {
    use super::*;

    #[test]
    fn test_uom_display() {
        let length = Length::new::<meter>(1.0);
        println!(
            "{}",
            length.into_format_args(meter, uom::fmt::DisplayStyle::Abbreviation)
        );

        println!(
            "{}",
            MegatonsPerYear::new(1.0)
                ._value
                .into_format_args(megaton_per_year, uom::fmt::DisplayStyle::Abbreviation)
        );
    }
}
