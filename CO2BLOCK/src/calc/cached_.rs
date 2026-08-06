use std::{
    collections::HashMap,
    sync::{LazyLock, RwLock},
};

use uom::si::f64::{Length, Ratio};

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

pub static CALC: LazyLock<Nordbotten> = LazyLock::new(Default::default);

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
        let result = super::nordbotten_impl(r, R, psi, R_ext, gamma);
        self.cache.write().unwrap().insert(args, result);
        result
    }
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
