use std::{marker::PhantomData, str::FromStr};
use uom::{
    Kind,
    si::{
        Dimension, Quantity, Units, amount_of_substance::mole, electric_current::ampere,
        length::meter, luminous_intensity::candela, mass::kilogram,
        thermodynamic_temperature::kelvin, time::second,
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

pub use mass_rate::megaton_per_year;
pub mod mass_rate {
    uom::unit! {
        system: uom::si;
        quantity: uom::si::mass_rate;
        @megaton_per_year: prefix!(mega) * 1.0_E3 / 3.1536_E7; "Mt/y", "megaton per year", "megatons per year";
    }
}

pub mod permeability_system {
    #[macro_use]
    mod length {
        uom::quantity! {
            quantity: Length; "length";
            dimension: Q<P1>;
            units {
                @meter: 1.0_E0; "m", "meter", "meters";
            }
        }
    }

    #[macro_use]
    pub mod permeability_quantity {
        uom::quantity! {
            quantity: Permeability; "permeability";
            dimension: Q<P2>;
            units {
                @square_meter: 1.0_E0; "m^2", "square meter", "square meters";
                @milli_darcy: prefix!(milli) * 9.86923_E-13; "mD", "millidarcy", "millidarcies";
            }
        }
    }

    uom::system! {
        quantities: Q {
            length: meter, L;
        }
        units: U {
            mod length::Length,
            mod permeability_quantity::Permeability,
        }
    }

    pub mod f64 {
        mod permeability {
            pub use super::super::*;
        }

        Q!(self::permeability, f64);
    }

    pub use f64::Permeability as PermeabilityQuantity;
    pub use permeability_quantity::milli_darcy;
    pub use permeability_quantity::square_meter as permeability_square_meter;
}

pub type Permeability = uom::si::f64::Area;

pub use permeability_system::{PermeabilityQuantity, milli_darcy, permeability_square_meter};

impl From<PermeabilityQuantity> for Permeability {
    fn from(value: PermeabilityQuantity) -> Self {
        Self {
            dimension: PhantomData,
            units: PhantomData,
            value: value.value,
        }
    }
}

impl From<Permeability> for PermeabilityQuantity {
    fn from(value: Permeability) -> Self {
        Self {
            dimension: PhantomData,
            units: PhantomData,
            value: value.value,
        }
    }
}

pub fn parse_permeability(
    value: &str,
) -> Result<Permeability, <PermeabilityQuantity as FromStr>::Err> {
    PermeabilityQuantity::from_str(value).map(Into::into)
}

pub fn format_permeability_mdarcy(value: Permeability) -> String {
    let permeability: PermeabilityQuantity = value.into();
    let fmt_args =
        PermeabilityQuantity::format_args(milli_darcy, uom::fmt::DisplayStyle::Abbreviation);
    format!("{}", fmt_args.with(permeability))
}
