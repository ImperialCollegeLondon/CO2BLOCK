#![deny(clippy::arithmetic_side_effects, clippy::as_conversions)]

#[macro_use]
extern crate uom;

mod calc;
mod mex;

// TODO: move MEX bindings to a separate crate
// pub use mex::nordbotten;
