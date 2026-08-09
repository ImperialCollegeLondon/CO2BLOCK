use uom::{
    si::{
        amount_of_substance::{self, mole},
        dynamic_viscosity::pascal_second,
        f64::*,
        mass::kilogram,
        mass_density::kilogram_per_cubic_meter,
        molar_volume::cubic_meter_per_mole,
        pressure::bar,
        ratio::ratio,
        thermodynamic_temperature::{degree_celsius, kelvin},
    },
    typenum::*,
};

pub fn brine_viscosity(temperature: ThermodynamicTemperature, salinity: Ratio) -> DynamicViscosity {
    let salinity_raw = salinity.get::<ratio>();
    let temperature_raw = temperature.get::<degree_celsius>();
    let brine_visc_raw = 1e-3
        * (0.1
            + 0.333 * salinity_raw
            + (1.65 + 91.9 * salinity_raw.powi(3))
                * (-(0.42 * (salinity_raw.powf(0.8) - 0.17).powi(2) + 0.045)
                    * temperature_raw.powf(0.8))
                .exp());

    DynamicViscosity::new::<uom::si::dynamic_viscosity::pascal_second>(brine_visc_raw)
}

pub fn gas_density(pressure: Pressure, temperature: ThermodynamicTemperature) -> MassDensity {
    let one_kelvin = ThermodynamicTemperature::new::<kelvin>(1.);
    let one_kelvin_dimless = one_kelvin / one_kelvin;
    let gas_constant = Energy::new::<uom::si::energy::joule>(8.314472)
        / one_kelvin_dimless
        / AmountOfSubstance::new::<amount_of_substance::mole>(1.);

    let p_c = Pressure::new::<bar>(73.773);
    let t_c = ThermodynamicTemperature::new::<kelvin>(304.128);

    let t_c_dimless = t_c / one_kelvin;

    let a =
        0.42748 * gas_constant.powi(P2::new()) * t_c_dimless.powi(P2::new()) * t_c_dimless.sqrt()
            / p_c;

    let b: MolarVolume = 0.08664 * gas_constant * t_c_dimless / p_c;

    let temp_dimless: Ratio = temperature / one_kelvin;

    let c2: MolarVolume = -gas_constant * temp_dimless / pressure;

    let c0 = -b * a / pressure / temp_dimless.sqrt();

    let c1 = b.powi(P2::new()) - b * c2 - a / pressure / temp_dimless.sqrt();

    let unit_molar_volume = MolarVolume::new::<cubic_meter_per_mole>(1.);

    // type check
    let _resid = unit_molar_volume.powi(P3::new())
        + unit_molar_volume.powi(P2::new()) * c2
        + unit_molar_volume * c1
        + c0;

    let c2: Ratio = c2 / unit_molar_volume;
    let c1: Ratio = c1 * unit_molar_volume.powi(N2::new());
    let c0: Ratio = c0 * unit_molar_volume.powi(N3::new());

    let roots =
        roots::find_roots_cubic_normalized(c2.get::<ratio>(), c1.get::<ratio>(), c0.get::<ratio>());

    let root: f64 = match roots {
        roots::Roots::No(_) => panic!(),
        roots::Roots::One(roots) => roots.into_iter().find(|el| el >= &0.).unwrap(),
        roots::Roots::Two(roots) => roots.into_iter().find(|el| el >= &0.).unwrap(),
        roots::Roots::Three(roots) => roots.into_iter().find(|el| el >= &0.).unwrap(),
        roots::Roots::Four(roots) => roots.into_iter().find(|el| el >= &0.).unwrap(),
    };

    let root = MolarVolume::new::<cubic_meter_per_mole>(root);

    Mass::new::<kilogram>(0.044) / root / AmountOfSubstance::new::<mole>(1.)
}

pub fn gas_viscosity(
    temperature: ThermodynamicTemperature,
    gas_density: MassDensity,
) -> DynamicViscosity {
    let reference_temperature = ThermodynamicTemperature::new::<kelvin>(304.);

    let reference_density = MassDensity::new::<kilogram_per_cubic_meter>(468.);

    let temp_normalized = temperature / reference_temperature;
    let dens_normalized = gas_density / reference_density;

    let unit = Ratio::new::<ratio>(1.);

    let mu0: Ratio = temp_normalized.sqrt()
        * (27.2246461 * unit - 16.6346068 / temp_normalized
            + 4.66920556 / (temp_normalized.powi(P2::new())))
        * 1e-6;
    let mu0 = DynamicViscosity::new::<pascal_second>(mu0.get::<ratio>());

    let exponent: Ratio = {
        let a10 = 0.248566120 * unit;
        let a11 = 0.004894942 * unit;
        let a20 = -0.373300660 * unit;
        let a21 = 1.22753488 * unit;
        let a30 = 0.363854523 * unit;
        let a31 = -0.774229021 * unit;
        let a40 = -0.0639070755 * unit;
        let a41 = 0.142507049 * unit;

        a10 * dens_normalized
            + a11 * dens_normalized / temp_normalized
            + a20 * dens_normalized.powi(P2::new())
            + a21 * dens_normalized.powi(P2::new()) / temp_normalized
            + a30 * dens_normalized.powi(P3::new())
            + a31 * dens_normalized.powi(P3::new()) / temp_normalized
            + a40 * dens_normalized.powi(P4::new())
            + a41 * dens_normalized.powi(P4::new()) / temp_normalized
    };
    mu0 * exponent.exp()
}
