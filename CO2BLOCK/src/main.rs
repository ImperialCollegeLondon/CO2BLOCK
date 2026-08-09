use co2block::calc::{InputReservoirParams, default_props};
use uom::si::{
    dynamic_viscosity::{DynamicViscosity, pascal_second},
    f64::{Length, MassDensity, Pressure, ThermodynamicTemperature},
    length::meter,
    mass_density::kilogram_per_cubic_meter,
    pressure::{bar, megapascal},
    thermodynamic_temperature::degree_celsius,
};

fn main() {
    let mut data = InputReservoirParams::default();
    data.rock_compress = Some(co2block::calc::default_props::rock_compressibility());
    data.pressure_top = Some(Pressure::new::<bar>(10.));
    data.co2_density = Some(MassDensity::new::<kilogram_per_cubic_meter>(1.));

    data.co2_viscosity = Some(DynamicViscosity::new::<pascal_second>(0.1));

    data.salinity = Some(default_props::salinity());

    data.princ_stress_ratio = Some(default_props::stress_ratio());

    data.temperature_center = Some(
        ThermodynamicTemperature::new::<degree_celsius>(15.)
            + default_props::temperature_grad() * Length::new::<meter>(1000.),
    );

    data.total_max_princ_stress = Some(Pressure::new::<megapascal>(10.));

    data.rock_friction_angle = Some(default_props::rock_friction_angle());

    let json = serde_json::to_string_pretty(&data).unwrap();
    println!("{json}");
    {
        let mut wtr = csv::Writer::from_path("input-template.csv").unwrap();
        wtr.serialize(data).unwrap();
    }
    let mut rdr = csv::Reader::from_path("input-template.csv").unwrap();
    let _data: InputReservoirParams = rdr.deserialize().into_iter().next().unwrap().unwrap();
}
