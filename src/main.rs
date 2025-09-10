/*
    rustbl: A compressible boundary layer analysis code.

    References:
     - "Hypersonic and High Temperature Gas Dyanmics", John D. Anderson

    @author: Nick Gibbons and Peter Jacobs
*/

#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]

use num_complex::Complex64;
use num_complex::ComplexFloat;
use std::ops::{Add,Mul,Div,Sub};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::io::prelude::*;
use std::env;

extern crate yaml_rust;
use yaml_rust::{YamlLoader};

pub mod state;
use crate::state::State;


// Super traits to allow mixed complex/real arithmetic
trait Cplx<T>: Mul<f64, Output = T> + Add<f64, Output = T>   + Div<f64, Output = T>   + Sub<f64, Output = T> {}
impl<T> Cplx<T> for T where T: Mul<f64, Output = T>   + Add<f64, Output = T>   + Div<f64, Output = T>   + Sub<f64, Output = T> {}


trait Mxd<T>:    Mul<T,Output = T>   + Add<T,   Output = T>   + Div<T,   Output = T>   + Sub<T,   Output = T>
               + Mul<f64, Output = f64> + Add<f64, Output = f64> + Div<f64, Output = f64> + Sub<f64, Output = f64> {}
impl<T> Mxd<T> for f64 where f64: Mul<T, Output = T>   + Add<T, Output = T>   + Div<T, Output = T>   + Sub<T, Output = T>
                                + Mul<f64, Output = f64> + Add<f64, Output = f64> + Div<f64, Output = f64> + Sub<f64, Output = f64> {}

#[derive(Debug)]
struct Config {
    R: f64,
    gamma: f64,
    Pr: f64,
    p_e: f64,
    u_e: f64,
    T_e: f64,
    T_wall: f64,
    x: f64,
}

fn read_config_file(filename: &str) -> Config {
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

fn rkf45_step(
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
    let k4 = f(t0 + 12.0*h/13.0, y0 + 1932.0*h*k1/2197.0 - 7200.0*h*k2/2197.0 + 7296.0*h*k3/2197.0, &pm);
    let k5 = f(t0 + h, y0 + 439.0*h*k1/216.0 - 8.0*h*k2 + 3680.0*h*k3/513.0 - 845.0*h*k4/4104.0, &pm);
    let k6 = f(t0 + h/2.0, y0 - 8.0*h*k1/27.0 + 2.0*h*k2 - 3544.0*h*k3/2565.0 + 1859.0*h*k4/4104.0 - 11.0*h*k5/40.0, &pm);
    // Now, do the integration as a weighting of the sampled data.
    let y1 = y0 + 16.0*h*k1/135.0 + 6656.0*h*k3/12825.0 + 28561.0*h*k4/56430.0 - 9.0*h*k5/50.0 + 2.0*h*k6/55.0;
    let err = (h*k1/360.0 - 128.0*h*k3/4275.0 - 2197.0*h*k4/75240.0 + h*k5/50.0 + 2.0*h*k6/55.0).abs();
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
    x: f64,
    T_wall: f64,
    h_wall: f64,
}

fn soft_max<T: ComplexFloat>(a: T, b: f64) -> T where T: Cplx<T>, f64: Mxd<T> {
    let da = a-b;
    let scale = 0.5*(a+b);
    let eps = 1e-3*scale + 1e-3;
    let max = scale + 0.5*ComplexFloat::sqrt(da*da + eps*eps);
    return max;
}

fn soft_max_derivative<T: ComplexFloat>(a: T, b: f64) -> T where T: Cplx<T>, f64: Mxd<T> {
    let da = a-b;
    let scale = 0.5*(a+b);
    let eps = 1e-3*scale + 1e-3;
    let sqrtval = ComplexFloat::sqrt(da*da + eps*eps);
    return 0.5 + 0.25/sqrtval*(2.0*da + 2.0*eps*1e-3/2.0);
}


const MU_REF: f64 = 1.716e-05;
const T_REF: f64 = 273.0;
const S: f64 = 111.0;

// T here is a generic type, not the temperature!
fn sutherland_mu<T: ComplexFloat>(TEMP: T) -> T where T: Cplx<T>, f64: Mxd<T> {
    return MU_REF*ComplexFloat::sqrt(TEMP/T_REF)*(TEMP/T_REF)*(T_REF + S)/(TEMP + S);
}

fn sutherland_mu_derivative<T: ComplexFloat>(TEMP: T) -> T where T: Cplx<T>, f64: Mxd<T> {
    return MU_REF*(T_REF+S)*ComplexFloat::sqrt(TEMP/T_REF)*(3.0*S+TEMP)/(2.0*T_REF*(S+TEMP)*(S+TEMP));
}

fn density_viscosity_product<T: ComplexFloat>(g: T, pm: &Parameters) -> T where T: Cplx<T>, f64: Mxd<T> {
/*
    Ratio of density x viscosity product at a given point in the boundary layer
*/
   let Temp = g*pm.h_e/pm.C_p;
   let SoftTemp = soft_max(Temp, 60.0);
   let rho = pm.p_e/(pm.R*SoftTemp);
   let mu = sutherland_mu(SoftTemp);
   return rho*mu/(pm.rho_e*pm.mu_e);
}

fn density_viscosity_product_derivative2<T: ComplexFloat>(g: T, pm: &Parameters) -> T where T: Cplx<T>, f64: Mxd<T> {
   let deltag = 0.0001 * g; // something not too big, not too small
   let C0 = density_viscosity_product(g, &pm);
   let C1 = density_viscosity_product(g+deltag, &pm);
   let dCdg = (C1-C0)/deltag;
   return dCdg;
}

fn density_viscosity_product_derivative<T: ComplexFloat>(g: T, pm: &Parameters) -> T where T: Cplx<T>, f64: Mxd<T> {
   let Temp = g*pm.h_e/pm.C_p;
   let SoftTemp = soft_max(Temp, 60.0);
   let rho = pm.p_e/(pm.R*SoftTemp);
   let mu = sutherland_mu(SoftTemp);

   let dTdg = pm.h_e/pm.C_p;
   let dSTdT = soft_max_derivative(Temp, 60.0);
   let dmudST = sutherland_mu_derivative(SoftTemp);
   let dmudT = dmudST*dSTdT;
   let dmudg  = dmudT*dTdg;

   let drhodST = -pm.p_e/pm.R/SoftTemp/SoftTemp;
   let drhodT = drhodST*dSTdT;
   let drhodg = drhodT*dTdg;

   return (rho*dmudg + mu*drhodg)/(pm.rho_e*pm.mu_e);
}

fn self_similar_ode(_t: f64, z: State, pm: &Parameters) -> State {
    let f = z.f; let fd = z.fd; let fdd = z.fdd;
    let g = z.g; let gd = z.gd; let y   = z.y;

    let C = density_viscosity_product(g, &pm);
    let dCdg= density_viscosity_product_derivative(g, &pm);
    let Cd = dCdg*gd; // oops it's dCdeta = dCdg*dgdeta

    let fddd = 1.0/C*(-f*fdd - Cd*fdd);
    let gdd = pm.Pr/C*(-gd*(Cd/pm.Pr+f) - C*pm.u_e*pm.u_e/pm.h_e*fdd*fdd);
    let yd = f64::sqrt(2.0*pm.xi)/pm.u_e*pm.h_e/pm.p_e*(pm.gamma-1.0)/pm.gamma*g;
    let dzdeta = State {f: fd, fd: fdd, fdd: fddd, g: gd, gd: gdd, y: yd};
    //println!("        Called ODE: g {:#} dCdg {:#} dCdg2 {:#}", g.re, dCdg.re, dCdg2.re);

    return dzdeta;
}

const NSTEPS: usize = 500;
fn integrate_through_bl(fdd: Complex64, gd: Complex64, pm: &Parameters) -> Vec<State> {
    let f  = Complex64::new(0.0,0.0);
    let fd = Complex64::new(0.0,0.0);
    let g = Complex64::new(pm.h_wall/pm.h_e, 0.0);
    let y = Complex64::new(0.0, 0.0);

    let mut eta0 = 0.0;
    let eta_final = 5.0;
    let deta = (eta_final-eta0)/(NSTEPS as f64);
    let mut z0 = State {f: f, fd: fd, fdd: fdd, g: g, gd: gd, y: y};
    let mut zs = Vec::with_capacity(NSTEPS+1);
    zs.push(z0.clone());

    for step in 0 .. NSTEPS {
        let (eta1, z1, _err) = rkf45_step(self_similar_ode, eta0, deta, z0, &pm);
        eta0 = eta1; z0 = z1;
        //let hi = z0.g.re*pm.h_e; let Ti = hi / pm.C_p;
        //let softTi = soft_max(Ti, 60.0);
        //println!("    step {:?}: [{:?},{:?},{:?},{:?},{:?},{:?}] T: {:?} ({:?})", step, z0.f.re, z0.fd.re, z0.fdd.re, z0.g.re, z0.gd.re, z0.y.re, Ti, softTi);
        //println!("           im: [{:?},{:?},{:?},{:?},{:?},{:?}] ", z0.f.im, z0.fd.im, z0.fdd.im, z0.g.im, z0.gd.im, z0.y.im);
        zs.push(z0.clone())
    }

    return zs;
}

fn skin_friction(z: State, pm: &Parameters) -> f64 {
/*
    Return tau, using equations 6.71 and 6.59 from Anderson
*/
    let rhomuw_on_rhomue = density_viscosity_product(z.g.re, &pm);
    let fddw = z.fdd.re;
    let Rex = pm.rho_e*pm.u_e/pm.mu_e*pm.x;
    let cf = f64::sqrt(2.0)*rhomuw_on_rhomue*fddw/f64::sqrt(Rex);
    let tau = 0.5*cf*pm.rho_e*pm.u_e*pm.u_e;
    return tau;
}

fn heat_transfer(z: State, pm: &Parameters) -> f64 {
/*
    Return q, using equations 6.79 and ??? from Anderson.
*/
    let rhomuw_on_rhomue = density_viscosity_product(z.g.re, &pm);
    let gdw = z.gd.re;
    let Rex = pm.rho_e*pm.u_e/pm.mu_e*pm.x;
    let qw = pm.u_e*pm.rho_e/f64::sqrt(2.0)/pm.Pr*rhomuw_on_rhomue*pm.h_e*gdw/f64::sqrt(Rex);
    return qw;
}


fn main() {
    println!("rustbl: A compressible boundary layer analysis code.");
    let mut config_file_name = "test.yaml";
    let args: Vec<String> = env::args().collect();
    if args.len()>1 { config_file_name = args[1].as_str(); }

    let config = read_config_file(config_file_name);
    let R     = config.R;
    let gamma = config.gamma;
    let Pr = config.Pr;

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
        x: x,
        k_e: k_e,
        T_wall: T_wall,
        h_wall: h_wall,
    };


    let mut error = 1e99;
    let tol = 1e-10;
    let mut iterations = 0;
    let mut fdd = 0.5;
    let mut gd  = 1.0;

    while error>tol {
        let eps = 1e-16;
        let fdd_pfdd= Complex64::new(fdd, eps);
        let gd_pfdd = Complex64::new(gd, 0.0);
        let states_dfdd= integrate_through_bl(fdd_pfdd, gd_pfdd, &pm);
        let state_dfdd= states_dfdd.last().unwrap();
        let d_fd_dfdd = (state_dfdd.fd.im)/eps; // d_fd_dfdd == df1_dtau
        let d_g_dfdd = (state_dfdd.g.im)/eps;   // d_g_dfdd  == df2_dtau

        let fdd_pgd = Complex64::new(fdd, 0.0);
        let gd_pgd  = Complex64::new(gd, eps);
        let states_dgd = integrate_through_bl(fdd_pgd, gd_pgd, &pm);
        let state_dgd = states_dgd.last().unwrap();
        let d_fd_dgd = (state_dgd.fd.im)/eps; // d_fd_dgd == df1_dq
        let d_g_dgd = (state_dgd.g.im)/eps;   // d_g_dgd  == df2_dq

        let fd_err = state_dgd.fd.re - 1.0; // f1 == fd_err
        let g_err  = state_dgd.g.re - 1.0;  // f2 == g_err
        error = f64::sqrt(fd_err*fd_err + g_err*g_err);

        // Two equation Newton's Method has a 2x2 jacobian that can be
        // inverted analytically. Do this here to get fdd and gd corrections.
        let diff_fdd = (fd_err*d_g_dgd/d_fd_dgd - g_err)
                      /(d_g_dfdd - d_g_dgd*d_fd_dfdd/d_fd_dgd);
        let diff_gd  = (fd_err*d_g_dfdd/d_fd_dfdd - g_err)
                      /(d_g_dgd - d_g_dfdd*d_fd_dgd/d_fd_dfdd);
        fdd += diff_fdd;
        gd  += diff_gd;

        iterations += 1;
        if iterations>100 { panic!("Too many iterations of newton solve"); }
    }

    println!("Got fdd {:#?} gd {:#?} in {:?} iters", fdd, gd, iterations);
    let state_initial = State::wall_state(fdd, gd, pm.h_wall, pm.h_e);
    let fdd_final = Complex64::new(fdd, 0.0);
    let gd_final = Complex64::new(gd, 0.0);
    let states = integrate_through_bl(fdd_final, gd_final, &pm);
    let state_final = states.last().unwrap();

    println!("Final state {:#?}", state_final);
    println!("Init state {:#?}", state_initial);
    let tauw = skin_friction(state_initial, &pm);
    println!("Skin Friction: {:#?} N/m2", tauw);
    let qw = heat_transfer(state_initial, &pm);
    println!("Heat Transfer : {:#?} W/m2", qw);

    let filename = config_file_name.replace(".yaml", ".dat");
    let file = File::create(filename.as_str()).expect("Unable to open for writing");
    let mut buf =  BufWriter::new(file);
    buf.write(b"# y vel T rho p\n").expect("Unable to write to file");

    for i in 0 .. NSTEPS+1 {
        let zi = states[i];
        let h = zi.g.re*pm.h_e; let T = h / pm.C_p; let rho = pm.p_e/(pm.R*T);
        let u = zi.fd.re*pm.u_e; let y = zi.y.re; let p = pm.p_e;

        write!(buf, "{:16.16e} {:16.16e} {:16.16e} {:16.16e} {:16.16e}\n", y, u, T, rho, p).expect("Unable to write line to file");
    }

    println!("Done.");
}
