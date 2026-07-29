use std::{
    collections::HashMap,
    sync::{LazyLock, RwLock},
};

use uom::si::{f64::Length, f64::Ratio};

#[derive(Clone)]
pub(crate) struct NordbottenArgs(
    pub(crate) Length,
    pub(crate) Length,
    pub(crate) Length,
    pub(crate) Length,
    pub(crate) Ratio,
);

impl NordbottenArgs {
    fn get(&self) -> (Length, Length, Length, Length, Ratio) {
        (self.0, self.1, self.2, self.3, self.4)
    }
}

pub(crate) static CALC: LazyLock<Nordbotten> = LazyLock::new(Default::default);

#[derive(Default)]
pub(crate) struct Nordbotten {
    cache: RwLock<HashMap<NordbottenArgs, Ratio>>,
}

impl PartialEq for NordbottenArgs {
    fn eq(&self, other: &Self) -> bool {
        self.0.value.to_bits() == other.0.value.to_bits()
            && self.1.value.to_bits() == other.1.value.to_bits()
            && self.2.value.to_bits() == other.2.value.to_bits()
            && self.3.value.to_bits() == other.3.value.to_bits()
            && self.4.value.to_bits() == other.4.value.to_bits()
    }
}
impl Eq for NordbottenArgs {}

impl std::hash::Hash for NordbottenArgs {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.value.to_bits().hash(state);
        self.1.value.to_bits().hash(state);
        self.2.value.to_bits().hash(state);
        self.3.value.to_bits().hash(state);
        self.4.value.to_bits().hash(state);
    }
}

impl Nordbotten {
    #[allow(non_snake_case)]
    pub(crate) fn eval(&self, args: NordbottenArgs) -> Ratio {
        if let Some(result) = self.cache.read().unwrap().get(&args) {
            return *result;
        }
        let (r, R, psi, R_ext, gamma) = args.get();
        let result = nordbotten_impl(r, R, psi, R_ext, gamma);
        self.cache.write().unwrap().insert(args, result);
        result
    }
}

#[allow(non_snake_case)]
fn nordbotten_impl(r: Length, R: Length, psi: Length, R_ext: Length, gamma: Ratio) -> Ratio {
    if r <= psi {
        return gamma * (psi / r).ln() + fd_nor(psi, R, R_ext);
    }

    if (r > psi) && (r <= R) {
        return fd_nor(r, R, R_ext);
    }

    Ratio::new::<uom::si::ratio::ratio>(0.)
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

#[test]
fn test_static() {
    use uom::si::length::meter;
    let args = NordbottenArgs(
        Length::new::<meter>(1.),
        Length::new::<meter>(2.),
        Length::new::<meter>(3.),
        Length::new::<meter>(4.),
        Ratio::new::<uom::si::ratio::ratio>(5.),
    );
    let calc = &*CALC;
    let result = calc.eval(args.clone());
    let result2 = calc.eval(args);
    assert_eq!(result, result2);
}
