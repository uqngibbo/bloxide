/*
    rustbl: A compressible boundary layer analysis code.

    @author: Nick Gibbons
*/

#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]

pub mod state;
use crate::state::State;

use ndarray::{array, Array1};

fn rkf45_step(
    f : fn(f64, Array1<f64>, &Parameters) -> Array1<f64>,
    t0 : f64,
    h : f64,
    y0 : Array1<f64>,
    pm: &Parameters
) -> (f64, Array1<f64>, Array1<f64>) {
    // Build up the sample point information as per the text book descriptions.
    let k1 = f(t0, y0.clone(), &pm);
    let k2 = f(t0 + h/4.0, y0.clone() + 0.25*h*k1.clone(), &pm);
    let k3 = f(t0 + 3.0*h/8.0, y0.clone() + 3.0*h*k1.clone()/32.0 + 9.0*h*k2.clone()/32.0, &pm);
    let k4 = f(t0 + 12.0*h/13.0, y0.clone() + 1932.0*h*k1.clone()/2197.0 - 7200.0*h*k2.clone()/2197.0 +
               7296.0*h*k3.clone()/2197.0, &pm);
    let k5 = f(t0 + h, y0.clone() + 439.0*h*k1.clone()/216.0 - 8.0*h*k2.clone() + 3680.0*h*k3.clone()/513.0 -
               845.0*h*k4.clone()/4104.0, &pm);
    let k6 = f(t0 + h/2.0, y0.clone() - 8.0*h*k1.clone()/27.0 + 2.0*h*k2.clone() - 3544.0*h*k3.clone()/2565.0 +
               1859.0*h*k4.clone()/4104.0 - 11.0*h*k5.clone()/40.0, &pm);
    // Now, do the integration as a weighting of the sampled data.
    let y1 = y0 + 16.0*h*k1.clone()/135.0 + 6656.0*h*k3.clone()/12825.0 +
        28561.0*h*k4.clone()/56430.0 - 9.0*h*k5.clone()/50.0 + 2.0*h*k6.clone()/55.0;
    let err = (h*k1/360.0 - 128.0*h*k3/4275.0 - 2197.0*h*k4/75240.0 +
               h*k5/50.0 + 2.0*h*k6/55.0).mapv(f64::abs);
    (t0+h, y1, err)
}

fn rkf45_step_2(
    f : fn(f64, State, &Parameters) -> State,
    t0 : f64,
    h : f64,
    y0 : State,
    pm: &Parameters
) -> (f64, State, State) {
    // Build up the sample point information as per the text book descriptions.
    let k1 = f(t0, y0, &pm);
    let k2 = f(t0 + h/4.0, y0 + 0.25*h*k1, &pm);
    let k3 = f(t0 + 3.0*h/8.0, y0 + 3.0*h*k1/32.0 + 9.0*h*k2/32.0, &pm);
    let k4 = f(t0 + 12.0*h/13.0, y0 + 1932.0*h*k1/2197.0 - 7200.0*h*k2/2197.0 +
               7296.0*h*k3/2197.0, &pm);
    let k5 = f(t0 + h, y0 + 439.0*h*k1/216.0 - 8.0*h*k2 + 3680.0*h*k3/513.0 -
               845.0*h*k4/4104.0, &pm);
    let k6 = f(t0 + h/2.0, y0 - 8.0*h*k1/27.0 + 2.0*h*k2 - 3544.0*h*k3/2565.0 +
               1859.0*h*k4/4104.0 - 11.0*h*k5/40.0, &pm);
    // Now, do the integration as a weighting of the sampled data.
    let y1 = y0 + 16.0*h*k1/135.0 + 6656.0*h*k3/12825.0 +
        28561.0*h*k4/56430.0 - 9.0*h*k5/50.0 + 2.0*h*k6/55.0;
    let err = (h*k1/360.0 - 128.0*h*k3/4275.0 - 2197.0*h*k4/75240.0 +
               h*k5/50.0 + 2.0*h*k6/55.0).abs();
    (t0+h, y1, err)
}


#[derive(Clone, Copy)]
struct Parameters {
    R: f64,
    gamma: f64,
    C_p: f64,
    Pr: f64,
    p_e: f64,
    T_e: f64,
    rho_e: f64,
    h_e: f64,
    mu_e: f64,
    u_e: f64,
    k_e: f64,
    xi: f64,
    T_wall: f64,
    h_wall: f64,
}

//use std::time::Instant;

const MU_REF: f64 = 1.716e-05;
const T_REF: f64 = 273.0;
const S: f64 = 111.0;

fn sutherland_mu(T: f64) -> f64 {
    return MU_REF*f64::sqrt(T/T_REF)*(T/T_REF)*(T_REF + S)/(T + S);
}

fn sutherland_mu_derivative(T: f64) -> f64 {
    return MU_REF*(T_REF+S)*f64::sqrt(T/T_REF)*(3.0*S+T)/(2.0*T_REF*(S+T)*(S+T));
}

fn density_viscosity_product(g: f64, pm: &Parameters) -> f64 {
/*
    Ratio of density x viscosity product at a given point in the boundary layer
*/
   let T = g*pm.h_e/pm.C_p; 
   //T = f64::max(T, 100.0); // Adds non analyticity FIXME????
   let rho = pm.p_e/(pm.R*T);
   let mu = sutherland_mu(T);
   return rho*mu/(pm.rho_e*pm.mu_e);
}

fn density_viscosity_product_derivative(g: f64, pm: &Parameters) -> f64{
   let T = g*pm.h_e/pm.C_p; 
   //T = f64::max(T, 100.0); // Adds non analyticity FIXME????
   let rho = pm.p_e/(pm.R*T);
   let mu = sutherland_mu(T);

   let dmudT = sutherland_mu_derivative(T);
   let drhodT = -pm.p_e/pm.R/T/T;
   let dTdg = pm.h_e/pm.C_p;

   let drhodg = drhodT*dTdg;
   let dmudg  = dmudT*dTdg;

   return (rho*dmudg + mu*drhodg)/(pm.rho_e*pm.mu_e);
}
    

fn main() {
    println!("rustbl: A compressible boundary layer analysis code.");

	let R     = 287.1;
	let gamma = 1.4;
    let Pr = 0.71;

	let p_e = 2.303e3;
	let u_e = 604.5;
	let T_e = 108.1;
    let T_wall = 307.0;

	let C_p = gamma/(gamma-1.0)*R;
    let h_e = C_p*T_e;
    let rho_e = p_e/(R*T_e);
    let mu_e = sutherland_mu(T_e);
    let k_e = mu_e*C_p/Pr;
    let h_wall = C_p*T_wall;

    let x = 2.1125;  // metres
    let xi = rho_e * u_e * mu_e * x;

    let pm = Parameters {
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
        k_e: k_e,
        T_wall: T_wall,
        h_wall: h_wall,
	};


    fn self_similar_ode(_t: f64, z: Array1<f64>, pm: &Parameters) -> Array1<f64> {
        let f = z[0]; let fd = z[1]; let fdd = z[2];
        let g = z[3]; let gd = z[4]; let y   = z[5];

        //println!("Called ode with input eta={:?} z={:?}", _t, z);
        let C = density_viscosity_product(g, &pm);
        let Cd= density_viscosity_product_derivative(g, &pm);
        //println!("    C: {:?} Cd: {:?}", C, Cd);

        let fddd = 1.0/C*(-f*fdd - Cd*fdd);
        let gdd = pm.Pr/C*(-gd*(Cd/pm.Pr+f) - C*pm.u_e*pm.u_e/pm.h_e*fdd*fdd);
        let yd = f64::sqrt(2.0*pm.xi)/pm.u_e*pm.h_e/pm.p_e*(pm.gamma-1.0)/pm.gamma*g;
        let dzdeta = array![fd, fdd, fddd, gd, gdd, yd];

        //println!("    Output: {:?}", dzdeta);
        return dzdeta;
    }

    fn self_similar_ode_2(_t: f64, z: State, pm: &Parameters) -> State {
        let f = z.f; let fd = z.fd; let fdd = z.fdd;
        let g = z.g; let gd = z.gd; let y   = z.y;

        //println!("Called ode with input eta={:?} z={:?}", _t, z);
        let C = density_viscosity_product(g, &pm);
        let Cd= density_viscosity_product_derivative(g, &pm);
        //println!("    C: {:?} Cd: {:?}", C, Cd);

        let fddd = 1.0/C*(-f*fdd - Cd*fdd);
        let gdd = pm.Pr/C*(-gd*(Cd/pm.Pr+f) - C*pm.u_e*pm.u_e/pm.h_e*fdd*fdd);
        let yd = f64::sqrt(2.0*pm.xi)/pm.u_e*pm.h_e/pm.p_e*(pm.gamma-1.0)/pm.gamma*g;
        let dzdeta = State {f: fd, fd: fdd, fdd: fddd, g: gd, gd: gdd, y: yd};

        //println!("    Output: {:?}", dzdeta);
        return dzdeta;
    }

    let f=0.0;
    let fd = 0.0;
    let g = pm.h_wall/pm.h_e;
    let y = 0.0;
    let fdd = 0.5;
    let gd = 1.0;

    let mut eta0 = 0.0;
    let mut eta0b = 0.0;
    let nsteps = 500;
    let eta_final = 5.0;
    let deta = (eta_final-eta0)/(nsteps as f64);
    let mut z0 = array![f, fd, fdd, g, gd, y];
    let mut z0b = State {f: f, fd: fd, fdd: fdd, g: g, gd: gd, y: y};


    for _ in 0 .. nsteps {
        let (eta1, z1, _err) = rkf45_step(self_similar_ode, eta0, deta, z0, &pm);
        eta0 = eta1; z0 = z1;
        let (eta1b, z1b, _errb) = rkf45_step_2(self_similar_ode_2, eta0b, deta, z0b, &pm);
        eta0b = eta1b; z0b = z1b;
    }

    println!("  final eta= {:?}", eta0);
    {
        let f = z0[0]; let fd = z0[1]; let fdd = z0[2];
        let g = z0[3]; let gd = z0[4]; let y   = z0[5];
        println!(" f= {:?} fd= {:?} fdd= {:?}", f, fd, fdd);
        println!(" g= {:?} gd= {:?}   y= {:?}", g, gd, y);
        println!("Error in final values of fd= {:?} g= {:?}", fd-1.0, g-1.0);
    }
    {
        let f = z0b.f; let fd = z0b.fd; let fdd = z0b.fdd;
        let g = z0b.g; let gd = z0b.gd; let y   = z0b.y;
        println!(" f= {:?} fd= {:?} fdd= {:?}", f, fd, fdd);
        println!(" g= {:?} gd= {:?}   y= {:?}", g, gd, y);
        println!("Error in final values of fd= {:?} g= {:?}", fd-1.0, g-1.0);
    }

    //let nstep = 10000;
    //let h = 1.0/(nstep as f64);
    //let mut x0 = array![0.0, 0.0, 0.0];
    //let start = Instant::now();
    //for _ in 0..nstep {
    //    let (t1, x1, _err) = rkf45_step(test_system1, t0, h, x0);
    //    t0 = t1; x0 = x1;
    //}
	//let elapsed = start.elapsed();
    //println!("  elapsed_time= {:?}", elapsed);
    //println!("  x           = {}", x0);
    //println!("  exact       = {:?}", solution1(t0));
    println!("Done.");
}
