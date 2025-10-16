/*
    Config file reading and associated data structures.

    @author: Nick Gibbons
*/

use std::fs::File;
use std::io::prelude::*;

use yaml_rust::{YamlLoader,Yaml};

#[derive(Debug)]
pub struct Config {
    pub R: f64,
    pub gamma: f64,
    pub Pr: f64,
    pub p_e: f64,
    pub u_e: f64,
    pub T_e: f64,
    pub T_wall: f64,
    pub x: f64,
}

fn coerce_to_f64(node: Yaml) -> f64 {
/*
    rust_yaml is very pedantic about type conversions and will not parse integer
    values using as_f64. This leads to trouble for users, who will hapilly enter
    values like "287" and expect them to be converted to 287.0 under the hood.

    This function uses a match statement to first attempt an f64 parse, then
    tries an i64 one that converts to f64 explicitly, then panics if neither
    operation is successful.
*/
    let val = node.as_f64();
    match val {
        Some(f64) => return val.unwrap(),
        None => return node.as_i64().expect("Failed to parse config file.") as f64,
    }
}

pub fn read_config_file(filename: &str) -> Config {
    let mut f = File::open(filename).expect(format!("Unable to open yaml file {}", filename).as_str());
    let mut buffer = String::new();
    f.read_to_string(&mut buffer).expect("Unable to parse file to string");

    let pages = YamlLoader::load_from_str(buffer.as_str()).unwrap();
    let cfg = &pages[0];

    let config = Config {
        R:      coerce_to_f64(cfg["R"].clone()),
        gamma:  coerce_to_f64(cfg["gamma"].clone()),
        Pr:     coerce_to_f64(cfg["Pr"].clone()),
        p_e:    coerce_to_f64(cfg["p_e"].clone()),
        u_e:    coerce_to_f64(cfg["u_e"].clone()),
        T_e:    coerce_to_f64(cfg["T_e"].clone()),
        T_wall: coerce_to_f64(cfg["T_wall"].clone()),
        x:      coerce_to_f64(cfg["x"].clone()),
    };
    return config
}

