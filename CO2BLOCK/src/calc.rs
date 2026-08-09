// #![allow(unused)]
mod args;
pub mod eos;
pub mod units;

pub use units::*;

pub use self::args::*;
use ::std::{
    f64::consts::{FRAC_PI_4, PI, TAU},
    ops::{Div, Mul},
};
use num_traits::{MulAdd, ToPrimitive, Zero};
use peroxide::{fuga::LambertWAccuracyMode, special::function::lambert_wm1};
use std::ops::Sub;
use uom::si::{Dimension, Quantity, Units, area::square_meter, f64::*, ratio::ratio};

// TODO: turn into a class with the fixed part computed in a constructor
pub fn co2block_with_placement(
    step: Length,
    num_steps: usize,
    reservoir: ReservoirParams,
    injection: InjectionParams,
    correction: Correction,
) -> MassRate {
    let num_wells = num_steps.mul_add(2, 1).pow(2);
    // WARNING: pay attention to where it's used
    let _guess_well_rate = {
        // WARNING: may not be a good guess
        let guess_well_rate_raw =
            reservoir.rock.permeability.get::<square_meter>() * 1e-13 / num_wells.to_f64().unwrap();
        let guess_well_mass_rate = MassRate::new::<megaton_per_year>(guess_well_rate_raw);
        guess_well_mass_rate / reservoir.gas.density
    };

    // NOTE: my new guess
    // WARNING: doen't this method require fixed-point iterations?
    let guess_well_rate = injection.max_well_rate / reservoir.gas.density / 10.;

    let nord_coef = NordbottenCoeff::from_args(NordbottenArgs {
        well_inj_rate: guess_well_rate,
        inj_time: injection.duration,
        porosity: reservoir.rock.porosity,
        thickness: reservoir.domain.thickness,
        visc_wat: reservoir.water.visc,
        visc_gas: reservoir.gas.visc,
        compr_rock: reservoir.rock.compress_rock,
        compr_wat: reservoir.water.compress,
        permeability: reservoir.rock.permeability,
        area: reservoir.domain.area,
    });

    let (_counts, coefs) = simulate_placement(step, num_steps, injection.well_radius, &nord_coef);

    let char_pressure = guess_well_rate * reservoir.water.visc
        / (reservoir.domain.thickness * reservoir.rock.permeability * TAU);

    let correction_term = if let Correction::On = correction {
        calc_correction(
            reservoir.water.visc,
            reservoir.gas.visc,
            num_wells as u64,
            calc_influence_radius(
                reservoir.rock.permeability,
                injection.duration,
                reservoir.water.visc,
                calc_total_compress(
                    reservoir.rock.compress_rock,
                    reservoir.rock.porosity,
                    reservoir.water.compress,
                ),
            ),
            calc_plume_ext(
                guess_well_rate,
                injection.duration,
                reservoir.rock.porosity,
                reservoir.domain.thickness,
            ),
            step,
        )
    } else {
        Zero::zero()
    };

    // all coefficients.sub(correction).mul(char_pressure)
    let over_pressure: Pressure = {
        let all_coefs: Ratio = coefs.into_iter().sum();
        all_coefs.sub(correction_term).mul(char_pressure)
    };

    let limit_rate = {
        let limit_pressure = calc_limit_pressure(
            reservoir.rock.max_principal_stress,
            reservoir.rock.top_pressure,
            reservoir.rock.friction_angle,
            reservoir.rock.stress_ratio,
            reservoir.rock.cohesion,
            reservoir.rock.tensile_strength,
        );
        calc_limit_rate(
            limit_pressure,
            over_pressure,
            guess_well_rate,
            reservoir.water.visc,
            reservoir.gas.visc,
            reservoir.rock.permeability,
            reservoir.domain.thickness,
            num_wells as u64,
            reservoir.gas.density,
        )
    };

    let updated: MassRate = update_rate(
        injection.max_well_rate,
        limit_rate,
        reservoir.gas.density,
        step,
        reservoir.rock.porosity,
        reservoir.domain.thickness,
        injection.duration,
        false,
    );

    updated
}

fn calc_theta(friction: Angle) -> Ratio {
    let fsin: Ratio = friction.sin();
    let one = Ratio::new::<ratio>(1.);
    (one - fsin) / (one + fsin)
}

fn calc_limit_pressure_shear(
    max_principal_stress: Pressure,
    top_pressure: Pressure,
    friction_angle: Angle,
    stress_ratio: Ratio,
    cohesion: Pressure,
) -> Pressure {
    let max_eff_princ_stress = max_principal_stress - top_pressure; //effective maximum principal stress [MPa]

    let theta = calc_theta(friction_angle);
    let one = Ratio::new::<ratio>(1.);
    let (fsin, fcos) = friction_angle.sin_cos();

    (stress_ratio - theta) / (one - theta) * max_eff_princ_stress + cohesion * fcos / fsin
}

fn calc_limit_pressure(
    max_principal_stress: Pressure,
    top_pressure: Pressure,
    friction_angle: Angle,
    stress_ratio: Ratio,
    cohesion: Pressure,
    tensile_strength: Pressure,
) -> Pressure {
    let shear_pres = calc_limit_pressure_shear(
        max_principal_stress,
        top_pressure,
        friction_angle,
        stress_ratio,
        cohesion,
    );

    let max_eff_princ_stress = max_principal_stress - top_pressure; //effective maximum principal stress [MPa]

    let s3 = stress_ratio * max_eff_princ_stress;

    let limit_pres_tensile = s3 + tensile_strength;

    shear_pres.min(limit_pres_tensile)
}

fn gamma(v_c: DynamicViscosity, v_w: DynamicViscosity) -> Ratio {
    v_c / v_w
}

fn delta(v_c: DynamicViscosity, v_w: DynamicViscosity) -> Ratio {
    Ratio::new::<ratio>(1.) - gamma(v_c, v_w)
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

fn _calc_max_num_wells(area: Area, dist_min: Length) -> uom::si::u32::Ratio {
    let res_float = area / (dist_min * dist_min);
    let res_raw: u32 = res_float
        .value
        .to_u32()
        .expect("divison + floor result must be valid & round");
    uom::si::u32::Ratio::new::<ratio>(res_raw)
}

fn _calc_well_dist_max(area: Area) -> Length {
    (area * 2.).sqrt() / 2.
}

// #[cached(
//     key = "(u64,u64,u64)",
//     convert = r#"{ calc_args_to_key(inter_well_dist, num_x, num_y) }"#
// )]
// TODO: parallelize
fn simulate_placement(
    step: Length,
    num_steps: usize,
    well_radius: Length,
    nord_coef: &NordbottenCoeff,
) -> (Vec<u64>, Vec<Ratio>) {
    let num_steps = {
        let max_num_steps = 2_f64
            .powf(-0.5)
            .mul(nord_coef.radius_reservoir / step)
            .floor::<ratio>()
            .get::<ratio>()
            .to_usize()
            .unwrap();
        num_steps.min(max_num_steps)
    };

    // compute all counts first
    let counts = {
        let mut counts = vec![0u64; num_steps + 1];
        counts[0] = 1;
        for num_steps_x in 1..=num_steps {
            for num_steps_y in 0..=num_steps_x {
                let count: u64 = match num_steps_y {
                    0 => 4,
                    x if (1..num_steps_x).contains(&x) => 8,
                    x if x == num_steps_x => 4,
                    _ => unreachable!("0 <= num_steps_y <= num_steps_x"),
                };
                counts[num_steps_x] += count;
            }
        }
        counts
    };

    let mut coefs = vec![Ratio::zero(); num_steps + 1];

    let c00 = nord_coef.calc(well_radius);
    coefs[0] += c00;

    for num_steps_x in 1..=num_steps {
        let x2 = num_steps_x * num_steps_x;
        for num_steps_y in 0..=num_steps_x {
            let y2 = num_steps_y * num_steps_y;
            let r2 = (x2 + y2)
                .to_f64()
                .expect("sum of integer squares should be representable in f64");

            let dist = step * r2.sqrt();

            let coef = nord_coef.calc(dist);

            let count: u64 = match num_steps_y {
                0 => 4,
                x if (1..num_steps_x).contains(&x) => 8,
                x if x == num_steps_x => 4,
                _ => unreachable!("0 <= num_steps_y <= num_steps_x"),
            };

            coefs[num_steps_x] += coef * count.to_f64().unwrap();
        }
    }
    (counts, coefs)
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

/*
let ans: Quantity<
    dyn Dimension<
        L = NInt<UInt<UInt<UInt<UTerm, B1>, B0>, B0>>,
        M = PInt<UInt<UTerm, B1>>,
        T = NInt<UInt<UTerm, B1>>,
        I = Z0,
        Th = Z0,
        N = Z0,
        J = Z0,
        Kind = dyn Kind + 'static
    >,
    dyn Units<
        f64,
        length = meter,
        mass = kilogram,
        time = second,
        electric_current = ampere,
        thermodynamic_temperature = kelvin,
        amount_of_substance = mole,
        luminous_intensity = candela
    >,
    f64
>
*/

fn reservoir_radius(area: Area) -> Length {
    (area / PI).sqrt()
}

pub struct NordbottenCoeff {
    radius_plume: Length,
    radius_influence: Length,
    radius_reservoir: Length,
    gas_to_water_visc: Ratio,
    big_influence_term: Ratio,
}

pub struct NordbottenArgs {
    well_inj_rate: VolumeRate,
    inj_time: Time,
    porosity: Ratio,
    thickness: Length,
    visc_wat: DynamicViscosity,
    visc_gas: DynamicViscosity,
    compr_rock: CompressibilityCoefficient,
    compr_wat: CompressibilityCoefficient,
    permeability: Area,
    area: Area,
}

impl NordbottenCoeff {
    fn from_args(args: NordbottenArgs) -> Self {
        let radius_plume = {
            let avg_plume_ext = calc_plume_ext(
                args.well_inj_rate,
                args.inj_time,
                args.porosity,
                args.thickness,
            );

            calc_equiv_plume_ext(avg_plume_ext, args.visc_wat, args.visc_gas)
        };

        let radius_influence = {
            let compr_total = args.compr_rock + args.porosity * args.compr_wat;
            calc_influence_radius(args.permeability, args.inj_time, args.visc_wat, compr_total)
        };

        Self::new(
            radius_plume,
            radius_influence,
            reservoir_radius(args.area),
            args.visc_gas / args.visc_wat,
        )
    }
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

fn calc_limit_rate(
    limit_pressure: Pressure,
    over_pressure: Pressure,
    guess_rate: VolumeRate,
    visc_w: DynamicViscosity,
    visc_g: DynamicViscosity,
    permeability: Area,
    res_thickness: Length,
    num_wells: u64,
    dens_g: MassDensity,
) -> MassRate {
    let b_term = (visc_w - visc_g) / (permeability * res_thickness)
        * (num_wells.to_f64().unwrap() / 4. + 1.)
        / (TAU * 2.);
    let exp_arg = -over_pressure / (guess_rate * b_term);
    let exp_mult = -limit_pressure / (guess_rate * b_term);
    let w: Ratio = exp_mult * exp_arg.exp();
    let lambert_sol = lambert_wm1(w.value, LambertWAccuracyMode::Precise);
    let limit_rate: VolumeRate = -limit_pressure / b_term / lambert_sol;
    limit_rate * dens_g
}

fn update_rate(
    max_rate: MassRate,
    limit_rate: MassRate,
    dens_g: MassDensity,
    step: Length,
    poro: Ratio,
    res_thickness: Length,
    inj_time: uom::si::f64::Time,
    is_central_well: bool,
) -> MassRate {
    let central_well_mult: f64 = if is_central_well { f64::INFINITY } else { 1. };
    let limit_rate_threshold: MassRate =
        central_well_mult * dens_g * step.powi(uom::typenum::P2::new()) * poro * res_thickness
            / inj_time
            * FRAC_PI_4;

    limit_rate.min(max_rate).min(limit_rate_threshold)
}

fn calc_plume_ext(
    well_inj_rate: VolumeRate,
    inj_time: uom::si::f64::Time,
    poro: Ratio,
    res_thickness: Length,
) -> Length {
    well_inj_rate
        .mul(inj_time)
        .div(poro * res_thickness * PI)
        .sqrt()
}

fn calc_omega(visc_w: DynamicViscosity, visc_g: DynamicViscosity) -> Ratio {
    (visc_g + visc_w) / (visc_g - visc_w) * (visc_g / visc_w).sqrt().ln() - Ratio::new::<ratio>(1.)
}

fn calc_equiv_plume_ext(
    plume_ext: Length,
    visc_w: DynamicViscosity,
    visc_g: DynamicViscosity,
) -> Length {
    calc_omega(visc_w, visc_g).exp() * plume_ext
}

struct HashedF64 {
    val: f64,
}

impl std::hash::Hash for HashedF64 {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.val.to_bits().hash(state);
    }
}

impl From<f64> for HashedF64 {
    fn from(val: f64) -> Self {
        Self { val }
    }
}

impl From<HashedF64> for f64 {
    fn from(value: HashedF64) -> Self {
        value.val
    }
}

impl<D: Dimension, U: Units<f64>> From<Quantity<D, U, f64>> for HashedF64 {
    fn from(val: Quantity<D, U, V>) -> Self {
        val.value.into()
    }
}

mod tests {
    #[test]
    fn test_uom_display() {
        use uom::si::{f64::Length, length::meter};
        let length = Length::new::<meter>(1.0);
        println!(
            "{}",
            length
                .into_format_args(meter, uom::fmt::DisplayStyle::Abbreviation)
                .to_string()
        );
    }
}
