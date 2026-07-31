// use crate::calc::*;
// use rustmex::{
//     Missing,
//     convert::{IntoRust, ToMatlab},
//     mxArray,
// };
// use uom::si::{f64::Length, length::meter};

// #[rustmex::entrypoint]
// #[allow(non_snake_case)]
// pub fn nordbotten(mut lhs: rustmex::LhsAns, _rhs: rustmex::Rhs) -> rustmex::Result<()> {
//     let r = parse_len(_rhs.first())?;
//     let R = parse_len(_rhs.get(1))?;
//     let psi = parse_len(_rhs.get(2))?;
//     let R_ext = parse_len(_rhs.get(3))?;
//     let gamma: uom::si::f64::Ratio = parse_raw(_rhs.get(4))?.into();
//     let result = nordbotten_impl(r, R, psi, R_ext, gamma);
//     *lhs.ans_mut() = Some(result.value.to_matlab());
//     Ok(())
// }
// fn parse_raw(input: Option<&&mxArray>) -> Result<f64, rustmex::Error> {
//     let val = input.error_if_missing("0", "missing variable")?.to_rust()?;
//     Ok(val)
// }

// fn parse_len(input: Option<&&mxArray>) -> Result<Length, rustmex::Error> {
//     let val = parse_raw(input)?;
//     let test: Length = Length::new::<meter>(val);
//     Ok(test)
// }
