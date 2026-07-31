use std::num::NonZeroU64;

use uom::si::{
    f64::{Length, Ratio},
    length::meter,
    mass_rate::ton_per_day,
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
    pub area: uom::si::f64::Area,
}

pub struct MegatonsPerYear {
    _value: uom::si::f64::MassRate,
}

impl MegatonsPerYear {
    pub fn new(value: f64) -> Self {
        let day_to_year = uom::si::f64::Time::new::<year>(1.).get::<day>();
        MegatonsPerYear {
            _value: uom::si::f64::MassRate::new::<ton_per_day>(value * prefix!(mega) / day_to_year),
        }
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

pub fn calculate(_args: Args) {
    let _ = nordbotten_impl;
    todo!();
}
