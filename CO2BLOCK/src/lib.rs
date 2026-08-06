#![deny(clippy::arithmetic_side_effects, clippy::as_conversions)]

#[macro_use]
extern crate uom;

mod calc;
mod mex;
// pub use mex::nordbotten;
pub use calc::calculate;
