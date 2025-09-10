/*
    Storage area for things that do not change over the solving the boundary layer problem.

    @author: Nick Gibbons
*/

use crate::config::{Config};
use crate::viscosity::{sutherland_mu};

#[derive(Clone, Copy)]
pub struct Parameters {
    pub R: f64,
    pub gamma: f64,
    pub C_p: f64,
    pub Pr: f64,
    pub p_e: f64,
    pub T_e: f64,
    pub rho_e: f64,
    pub h_e: f64,
    pub mu_e: f64,
    pub u_e: f64,
    pub k_e: f64,
    pub xi: f64,
    pub x: f64,
    pub T_wall: f64,
    pub h_wall: f64,
}

impl Parameters{
    pub fn new(config: &Config) -> Self {
        let R     = config.R;
        let gamma = config.gamma;
        let Pr    = config.Pr;

        let p_e = config.p_e;
        let u_e = config.u_e;
        let T_e = config.T_e;
        let T_wall = config.T_wall;
        let x = config.x;  // metres

        let C_p = gamma/(gamma-1.0)*R;
        let h_e = C_p*T_e;
        let rho_e = p_e/(R*T_e);
        let mu_e = sutherland_mu(T_e);

        let k_e = mu_e*C_p/Pr;
        let h_wall = C_p*T_wall;

        let xi = rho_e * u_e * mu_e * x;
        Self {
            R: R,
            gamma: gamma,
            Pr: Pr,
            C_p: C_p,
            p_e: p_e,
            T_e: T_e,
            rho_e: rho_e,
            h_e: h_e,
            mu_e: mu_e,
            u_e: u_e,
            xi: xi,
            x: x,
            k_e: k_e,
            T_wall: T_wall,
            h_wall: h_wall,
        }
    }
}
