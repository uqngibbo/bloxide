/*
    Config file reading and associated data structures.

    @author: Nick Gibbons
*/

use std::fs::File;
use std::io::prelude::*;

use yaml_rust::{YamlLoader};

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

pub fn read_config_file(filename: &str) -> Config {
    let mut f = File::open(filename).expect(format!("Unable to open yaml file {}", filename).as_str());
    let mut buffer = String::new();
    f.read_to_string(&mut buffer).expect("Unable to parse file to string");

    let pages = YamlLoader::load_from_str(buffer.as_str()).unwrap();
    let cfg = &pages[0];

    let config = Config {
        R:      cfg["R"].as_f64().unwrap(),
        gamma:  cfg["gamma"].as_f64().unwrap(),
        Pr:     cfg["Pr"].as_f64().unwrap(),
        p_e:    cfg["p_e"].as_f64().unwrap(),
        u_e:    cfg["u_e"].as_f64().unwrap(),
        T_e:    cfg["T_e"].as_f64().unwrap(),
        T_wall: cfg["T_wall"].as_f64().unwrap(),
        x:      cfg["x"].as_f64().unwrap(),
    };
    return config
}

