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
    pub mod permeability_quantity {
        uom::quantity! {
            quantity: Permeability; "permeability";
            dimension: Q<P1>;
            units {
                @square_meter: 1.0_E0; "m^2", "square meter", "square meters";
                @milli_darcy: prefix!(milli) * 9.86923_E-13; "mD", "millidarcy", "millidarcies";
            }
        }
    }

    uom::system! {
        quantities: Q {
            permeability_quantity: square_meter, P;
        }
        units: U {
            mod permeability_quantity::Permeability,
        }
    }

    pub mod f64 {
        mod permeability {
            pub use super::super::*;
        }

        Q!(self::permeability, f64);
    }

    pub use f64::Permeability;
    pub use permeability_quantity::milli_darcy;
    pub use permeability_quantity::square_meter as permeability_square_meter;
}

pub use permeability_system::{Permeability, milli_darcy, permeability_square_meter};
