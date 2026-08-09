#![deny(clippy::arithmetic_side_effects, clippy::as_conversions)]

pub mod calc;

// TODO: move MEX bindings to a separate crate
mod mex;

#[macro_use]
extern crate uom;

mod pressure_gradient {
    use uom::{
        Kind,
        si::{
            Dimension, Quantity, Units, amount_of_substance::mole, electric_current::ampere,
            length::meter, luminous_intensity::candela, mass::kilogram,
            temperature_interval::kelvin, time::second,
        },
        typenum::*,
    };

    pub type PressureGradient = Quantity<
        dyn Dimension<
                L = N2,
                M = P1,
                T = N2,
                I = Z0,
                Th = Z0,
                N = Z0,
                J = Z0,
                Kind = dyn Kind + 'static,
            >,
        dyn Units<
                f64,
                length = meter,
                mass = kilogram,
                time = second,
                electric_current = ampere,
                thermodynamic_temperature = kelvin,
                amount_of_substance = mole,
                luminous_intensity = candela,
            >,
        f64,
    >;
}
