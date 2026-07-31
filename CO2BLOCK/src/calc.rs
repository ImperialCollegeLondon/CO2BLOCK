use uom::si::{f64::Length, f64::Ratio};

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

struct Args {
    correction: Correction,
    inter_well_dist_min: Len<kilometer>,
    inter_well_dist_max: Option<Len<kilometer>>,
}

struct Len<Unit> {
    value: Length,
    _unit: std::marker::PhantomData<Unit>,
}
impl<Unit: uom::si::length::Unit + uom::Conversion<f64, T = f64>> Len<Unit> {
    fn new(value: f64) -> Self {
        Len {
            value: Length::new::<Unit>(value),
            _unit: std::marker::PhantomData,
        }
    }

    fn get(self) -> Length {
        self.value
    }
}

use uom::si::length::kilometer;

enum Correction {
    Off,
    On,
}

fn calculate(args: Args) {}
