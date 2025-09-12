/*
    bloxide: A compressible boundary layer analysis code.

    References:
     - "Hypersonic and High Temperature Gas Dyanmics", John D. Anderson

    @author: Nick Gibbons and Peter Jacobs
*/

#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]

use std::env;

use bloxide::config::{read_config_file};
use bloxide::parameters::{Parameters};
use bloxide::{solve_boundary_layer, skin_friction, heat_transfer, boundary_layer_size, write_dat_file, reynolds_number, solve_adiabatic_boundary_layer};

fn main() {
    println!("bloxide: A compressible boundary layer analysis code.");
    let mut config_file_name = "test.yaml";
    let args: Vec<String> = env::args().collect();
    if args.len()>1 { config_file_name = args[1].as_str(); }

    let config = read_config_file(config_file_name);
    let pm = Parameters::new(&config);
    println!("{:#?}", config);

    let states = solve_boundary_layer(&pm);
    let state_initial = states[0];
    let state_final = states.last().unwrap();

    let tauw = skin_friction(state_initial, &pm);
    let qw = heat_transfer(state_initial, &pm);
    let ybl = boundary_layer_size(&states).expect("Cannot find bl size");
    let Rex = reynolds_number(pm.rho_e, pm.u_e, pm.mu_e, pm.x);
    let Ret = reynolds_number(pm.rho_e, pm.u_e, pm.mu_e, ybl);

    println!("Skin Friction:  {:5.5} N/m2", tauw);
    println!("Heat Transfer : {:5.5} W/cm2", qw/1e4);
    println!("99.9% BL size : {:5.5} mm", ybl*1000.0);
    println!("Rex: {:5.5} million   Ret: {:5.5}", Rex/1e6, Ret);

    let adiabatic_states = solve_adiabatic_boundary_layer(&pm);
    let adiabatic_state_initial = adiabatic_states[0];
    let hwall = adiabatic_state_initial.g.re*pm.h_e;
    let Twall = hwall / pm.C_p;
    println!("Adiabatic Wall Temp: {:5.5} K", Twall);

    let filename = config_file_name.replace(".yaml", ".dat");
    write_dat_file(states, filename.as_str(), &pm);

    println!("Done.");
}
