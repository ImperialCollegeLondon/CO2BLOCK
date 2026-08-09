#![warn(
    // clippy::arithmetic_side_effects,
    clippy::as_conversions
)]

pub mod calc;

// TODO: move MEX bindings to a separate crate
mod mex;

#[macro_use]
extern crate uom;
