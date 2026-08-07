#![allow(unused)]

use std::{
    num::NonZeroU64,
    ops::{Div, Mul},
};

use num_traits::{ToPrimitive, Zero};
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

mod cached_;

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

use cached::cached;

fn calc_args_to_key(
    inter_well_dist: Length,
    num_x: NonZeroU64,
    num_y: NonZeroU64,
) -> (u64, u64, u64) {
    (inter_well_dist.value.to_bits(), num_x.into(), num_y.into())
}

// #[cached(
//     key = "(u64,u64,u64)",
//     convert = r#"{ calc_args_to_key(inter_well_dist, num_x, num_y) }"#
// )]
pub fn calculate(
    nord_coef: &NordbottenCoeff,
    step: Length,
    well_radius: Length,
    mass_rate: MassRate,
    visc_w: DynamicViscosity,
    res_thickness: Length,
    permeability: Area,
    gas_density: MassDensity,
) {
    let (counts, coefs) = {
        let max_num_steps = 2_f64
            .powf(-0.5)
            .mul(nord_coef.radius_reservoir / step)
            .floor::<ratio>()
            .get::<ratio>()
            .to_usize()
            .unwrap();

        // let's store partial sums
        let mut coefs = vec![Ratio::zero(); max_num_steps];
        let mut counts = vec![0u64; max_num_steps];

        let c00 = nord_coef.calc(well_radius);
        coefs[0] += c00;
        counts[0] += 1;

        for num_steps_x in 1..=max_num_steps {
            // add all previous to the new sum
            {
                let prev = coefs[num_steps_x - 1];
                coefs[num_steps_x] += prev;
            }
            {
                let prev = counts[num_steps_x - 1];
                counts[num_steps_x] += prev;
            }
            let x2 = num_steps_x * num_steps_x;
            for num_steps_y in 0..=num_steps_x {
                let y2 = num_steps_y * num_steps_y;
                let r2 = (x2 + y2)
                    .to_f64()
                    .expect("sum of integer squares should be representable in f64");

                let dist = step * r2.sqrt();

                let coef = nord_coef.calc(dist);

                // symmetry
                let count: u64 = match num_steps_y {
                    0 => 4,
                    x if (1..num_steps_x).contains(&x) => 8,
                    x if x == num_steps_x => 4,
                    _ => unreachable!("0 <= num_steps_y <= num_steps_x"),
                };

                coefs[num_steps_x] += coef * count.to_f64().unwrap();
                counts[num_steps_x] += count;
            }
        }
        (counts, coefs)
    };

    let volume_rate: VolumeRate = mass_rate / gas_density;
    let char_pressure =
        volume_rate * visc_w / (res_thickness * permeability * std::f64::consts::TAU);
}

fn calc_correction(
    visc_w: DynamicViscosity,
    visc_g: DynamicViscosity,
    num_wells: u64,
    radius_influence: Length,
    avg_plume_ext: Length,
    grid_step: Length,
) -> Ratio {
    radius_influence
        .mul(avg_plume_ext)
        .div(grid_step.powi(uom::typenum::P2::new()))
        .ln()
        .mul(delta(visc_g, visc_w))
        .mul(num_wells.to_f64().expect("should be representable in f64"))
        .div(4.)
}

pub struct NordbottenCoeff {
    radius_plume: Length,
    radius_influence: Length,
    radius_reservoir: Length,
    gas_to_water_visc: Ratio,
    big_influence_term: Ratio,
}

impl NordbottenCoeff {
    fn new(
        radius_plume: Length,
        radius_influence: Length,
        radius_reservoir: Length,
        gas_to_water_visc: Ratio,
    ) -> Self {
        let big_influence_term: Ratio = (radius_influence > radius_reservoir)
            .then(|| {
                let coef = radius_influence / radius_reservoir;
                (coef * coef).mul(8. / 9.) - Ratio::new::<ratio>(3. / 4.)
            })
            .unwrap_or_else(Zero::zero);
        Self {
            radius_plume,
            radius_influence,
            radius_reservoir,
            gas_to_water_visc,
            big_influence_term,
        }
    }

    fn radius_reservoir(&self) -> Length {
        self.radius_reservoir
    }

    fn calc(&self, distance: Length) -> Ratio {
        let inv_dist_normalized = self.radius_plume / distance.min(self.radius_plume);

        let plume_term = inv_dist_normalized.ln().mul(self.gas_to_water_visc);

        let influence_term = self.fd_nor(distance.max(self.radius_plume));

        plume_term + influence_term
    }

    fn fd_nor(&self, distance: Length) -> Ratio {
        let limit_dist = self.radius_influence.min(self.radius_reservoir);

        let inv_dist_norm = limit_dist / distance.max(self.radius_reservoir);

        let dist_term = inv_dist_norm.ln();

        dist_term + self.big_influence_term * f64::from(distance < self.radius_reservoir)
    }
}

use peroxide::special::function::lambert_wm1;

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
