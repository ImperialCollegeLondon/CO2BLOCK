use uom::{
    Kind,
    si::{
        Dimension, Quantity, Units, amount_of_substance::mole, electric_current::ampere,
        length::meter, luminous_intensity::candela, mass::kilogram, temperature_interval::kelvin,
        time::second,
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

pub use permeability::milli_darcy;

pub mod permeability {
    uom::unit! {
        system: uom::si;
        quantity: uom::si::area;
        @milli_darcy: prefix!(milli) * 9.86923 * 10e-13; "mD", "millidarcy", "millidarcies";
    }
}
