#![doc = r#"
Original function

```matlab
% This function returns brine viscosity and CO2 density and viscosity.
% Brine viscosity is calculated according to Batzle and Wang (1992)
% CO2 density is calculated according to Redlich and  Kwong (1949) (with  the
% parameters proposed by Spycher et al. (2003).
% CO2 viscosity is calculated according to Altunin and Sakhabetdinov (1972)

function [brineviscosity, co2density, co2viscosity] = eos(T,p,salinity,co2dens)   % T in ºC, p in MPa, salinity in [ppm/1e6]

    % brine viscosity
    brineviscosity =(0.1+0.333*salinity+(1.65+91.9*salinity^3)...
        *exp(-(0.42*(salinity^0.8-0.17)^2+0.045)*T^0.8))/1e3;  %[Pa.s]


    T = T+273.15;               % transform [C] to [K]
    p = p*1e6 ;                 % transform [MPa] to [Pa]

    % CO2 density
    syms x
    a0 = 7.54;                  % constant [Pa m6 K^0.5 mol^-2]
    a1 = -4.13*10^-3;           % constant [Pa m6 K^0.5 mol^-2]
    b = 2.78*10^-5 ;            % constant [m3/mol]
    a = a0+a1*T ;
    R = 8.314472 ;              % gas constant [m3 Pa K?1 mol?1]
    eqn = x.^3-(R*T/p)*x.^2-(R*T*b/p-a/p/sqrt(T)+b.^2)*x - (a*b/p/sqrt(T)) == 0 ;
    V = vpa(solve(eqn,x,'real',true)) ;      % molar volume [m3/mol]
    co2density = double(0.044 / V );         % [kg/m3]

    % CO2 viscosity
    a10 = 0.248566120;
    a11 = 0.004894942 ;
    a20 = -0.373300660;
    a21 = 1.22753488 ;
    a30 = 0.363854523 ;
    a31 = -0.774229021 ;
    a40 = -0.0639070755 ;
    a41 = 0.142507049 ;
    Tr = T/304;                                    % reduced temperature
    dens_r = co2dens/468  ;                     % reduced density
    mu_0 = Tr^0.5* (27.2246461 - 16.6346068 /Tr + 4.66920556/(Tr^2))*1e-6 ; %[Pa s]
    co2viscosity = double(mu_0*exp(a10*dens_r + a11*dens_r/Tr + a20*dens_r^2 + a21*dens_r^2/Tr + ...
                a30*dens_r^3 + a31*dens_r^3/Tr  +  a40*dens_r^4 + a41*dens_r^4/Tr ));   %[Pa s]

end

```
"#]

use uom::{
    si::{
        amount_of_substance::{self, mole},
        f64::{Length, *},
        length::meter,
        mass::kilogram,
        mass_density::kilogram_per_cubic_meter,
        molar_volume::cubic_meter_per_mole,
        pressure::pascal,
        ratio::ratio,
        temperature_interval::terakelvin,
        thermodynamic_temperature::{degree_celsius, kelvin},
    },
    typenum::*,
};

struct ThermodynmacState {
    brine_viscosity: DynamicViscosity,
    gas_density: MassDensity,
    gas_viscosity: DynamicViscosity,
}

impl ThermodynmacState {
    fn from_eos(
        pressure: Pressure,
        temperature: ThermodynamicTemperature,
        salinity: Ratio,
    ) -> ThermodynmacState {
        let salinity_raw = salinity.get::<ratio>();
        let temperature_raw = temperature.get::<degree_celsius>();
        let brine_visc_raw = 1e-3
            * (0.1
                + 0.333 * salinity_raw
                + (1.65 + 91.9 * salinity_raw.powi(3))
                    * (-(0.42 * (salinity_raw.powf(0.8) - 0.17).powi(2) + 0.045)
                        * temperature_raw.powf(0.8))
                    .exp());

        let brine_viscosity =
            DynamicViscosity::new::<uom::si::dynamic_viscosity::pascal_second>(brine_visc_raw);

        let gas_constant = 8.314472; // J / (mol * K)
        let gas_constant = Energy::new::<uom::si::energy::joule>(8.314472)
            / ThermodynamicTemperature::new::<kelvin>(1.)
            / AmountOfSubstance::new::<amount_of_substance::mole>(1.);

        let c3 = MolarVolume::new::<cubic_meter_per_mole>(1.).powi(P3::new());
        let c2: MolarVolume = gas_constant * temperature / pressure;

        let one_kelvin = ThermodynamicTemperature::new::<kelvin>(1.);
        let one_kelvin_dimless = one_kelvin / one_kelvin;
        let temp_dimless: Ratio = temperature / one_kelvin;

        // BUG
        let a0 = Pressure::new::<pascal>(7.54)
            * Length::new::<meter>(1.).powi(P6::new())
            * AmountOfSubstance::new::<mole>(1.).powi(N2::new())
            * one_kelvin_dimless.sqrt();

        let a1 = Pressure::new::<pascal>(-4.13e-3)
            * Length::new::<meter>(1.).powi(P6::new())
            * AmountOfSubstance::new::<mole>(1.).powi(N2::new())
            / one_kelvin_dimless.sqrt();

        let a = a0 / temp_dimless.sqrt() + a1 * temp_dimless.sqrt();
        let b = MolarVolume::new::<cubic_meter_per_mole>(2.7e-5);

        let c0 /* : MolarVolume^3 */ = b * a / pressure;

        let c1 = c2 * b + b.powi(P2::new()) - a / pressure;

        let unit_molar_volume = MolarVolume::new::<cubic_meter_per_mole>(1.);

        // type check
        let _resid = unit_molar_volume.powi(P3::new())
            + unit_molar_volume.powi(P2::new()) * c2
            + unit_molar_volume * c1
            + c0;

        let c2: Ratio = c2 / unit_molar_volume;
        let c1: Ratio = c1 * unit_molar_volume.powi(N2::new());
        let c0: Ratio = c0 * unit_molar_volume.powi(N3::new());

        let roots = roots::find_roots_cubic_normalized(
            c2.get::<ratio>(),
            c1.get::<ratio>(),
            c0.get::<ratio>(),
        );

        let root: f64 = match roots {
            roots::Roots::No(_) => panic!(),
            roots::Roots::One(roots) => roots.into_iter().filter(|el| el >= &0.).next().unwrap(),
            roots::Roots::Two(roots) => roots.into_iter().filter(|el| el >= &0.).next().unwrap(),
            roots::Roots::Three(roots) => roots.into_iter().filter(|el| el >= &0.).next().unwrap(),
            roots::Roots::Four(roots) => roots.into_iter().filter(|el| el >= &0.).next().unwrap(),
        };

        let root = MolarVolume::new::<cubic_meter_per_mole>(root);
        let gas_density = Mass::new::<kilogram>(0.044) / root / AmountOfSubstance::new::<mole>(1.);

        let gas_viscosity: DynamicViscosity = {
            let reference_temperature = ThermodynamicTemperature::new::<kelvin>(304.);

            let reference_density = MassDensity::new::<kilogram_per_cubic_meter>(468.);

            let temp_normalized = temperature / reference_temperature;
            let dens_normalized = gas_density / reference_density;

            let unit = Ratio::new::<ratio>(1.);

            let reference_viscosity: Ratio = temp_normalized.sqrt()
                * (27.2246461 * unit - 16.6346068 / temp_normalized
                    + 4.66920556 / (temp_normalized.powi(P2::new())))
                * 1e-6;
            todo!()
        };

        ThermodynmacState {
            brine_viscosity,
            gas_density,
            gas_viscosity: (),
        }
    }
}
