#![deny(clippy::arithmetic_side_effects, clippy::as_conversions)]

pub mod calc;

// TODO: move MEX bindings to a separate crate
mod mex;

#[macro_use]
extern crate uom;

#[macro_use]
mod permeability {
    uom::quantity! {
        quantity: Permeability; "permeability";
        dimension: Q<P2, Z0, Z0>;
        units {
            @milli_darcy: prefix!(milli) * 9.86923 * 10e-13; "mD", "millidarcy", "millidarcies";
        }
    }
}
