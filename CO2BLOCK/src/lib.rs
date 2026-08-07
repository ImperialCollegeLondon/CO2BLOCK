#![deny(clippy::arithmetic_side_effects, clippy::as_conversions)]

#[macro_use]
extern crate uom;

pub mod calc;

// TODO: move MEX bindings to a separate crate
mod mex;
