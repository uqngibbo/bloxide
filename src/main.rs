/*
    rustbl: A compressible boundary layer analysis code.

    @author: Nick Gibbons
*/

#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused_variables)]

use num_complex::Complex64;
use num_complex::ComplexFloat;
use std::ops::{Add,Mul,Div,Sub};


pub mod state;
use crate::state::State;


// Super traits to allow mixed complex/real arithmetic
trait Cplx<T>: Mul<f64, Output = T> + Add<f64, Output = T>   + Div<f64, Output = T>   + Sub<f64, Output = T> {}
impl<T> Cplx<T> for T where T: Mul<f64, Output = T>   + Add<f64, Output = T>   + Div<f64, Output = T>   + Sub<f64, Output = T> {}


trait Mxd<T>:    Mul<T,Output = T>   + Add<T,   Output = T>   + Div<T,   Output = T>   + Sub<T,   Output = T>
               + Mul<f64, Output = f64> + Add<f64, Output = f64> + Div<f64, Output = f64> + Sub<f64, Output = f64> {}
impl<T> Mxd<T> for f64 where f64: Mul<T, Output = T>   + Add<T, Output = T>   + Div<T, Output = T>   + Sub<T, Output = T>
                                + Mul<f64, Output = f64> + Add<f64, Output = f64> + Div<f64, Output = f64> + Sub<f64, Output = f64> {}


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
   let T = g*pm.h_e/pm.C_p; 
   //T = f64::max(T, 100.0); // Adds non analyticity FIXME????
   let rho = pm.p_e/(pm.R*T);
   let mu = sutherland_mu(T);
   return rho*mu/(pm.rho_e*pm.mu_e);
}

fn density_viscosity_product_derivative<T: ComplexFloat>(g: T, pm: &Parameters) -> T where T: Cplx<T>, f64: Mxd<T> {
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

    return dzdeta;
}

fn integrate_through_bl(fdd: Complex64, gd: Complex64, pm: &Parameters) -> State {
    let f  = Complex64::new(0.0,0.0);
    let fd = Complex64::new(0.0,0.0);
    let g = Complex64::new(pm.h_wall/pm.h_e, 0.0);
    let y = Complex64::new(0.0, 0.0);

    let mut eta0 = 0.0;
    let nsteps = 500;
    let eta_final = 5.0;
    let deta = (eta_final-eta0)/(nsteps as f64);
    let mut z0 = State {f: f, fd: fd, fdd: fdd, g: g, gd: gd, y: y};

    for step in 0 .. nsteps {
        let (eta1, z1, _err) = rkf45_step(self_similar_ode, eta0, deta, z0, &pm);
        eta0 = eta1; z0 = z1;
    }

    return z0;
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


    let mut error = 1e99;
    let tol = 1e-6;
    let mut iterations = 0;
    let mut fdd = 0.5;
    let mut gd  = 1.0;

    while error>tol {
        let eps = 1e-16;
        let fdd_pfdd= Complex64::new(fdd, eps);
        let gd_pfdd = Complex64::new(gd, 0.0);
        let state_dfdd= integrate_through_bl(fdd_pfdd, gd_pfdd, &pm);
        let d_fd_dfdd = (state_dfdd.fd.im)/eps; // d_fd_dfdd == df1_dtau
        let d_g_dfdd = (state_dfdd.g.im)/eps;   // d_g_dfdd  == df2_dtau

        let fdd_pgd = Complex64::new(fdd, 0.0);
        let gd_pgd  = Complex64::new(gd, eps);
        let state_dgd = integrate_through_bl(fdd_pgd, gd_pgd, &pm);
        let d_fd_dgd = (state_dgd.fd.im)/1e-16; // d_fd_dgd == df1_dq
        let d_g_dgd = (state_dgd.g.im)/1e-16;   // d_g_dgd  == df2_dq

        let fd_err = state_dgd.fd.re - 1.0; // f1 == fd_err
        let g_err  = state_dgd.g.re - 1.0;  // f2 == g_err
        error = f64::sqrt(fd_err*fd_err + g_err*g_err);

        let diff_fdd = (fd_err*d_g_dgd/d_fd_dgd - g_err)
                      /(d_g_dfdd - d_g_dgd*d_fd_dfdd/d_fd_dgd);
        let diff_gd  = (fd_err*d_g_dfdd/d_fd_dfdd - g_err)
                      /(d_g_dgd - d_g_dfdd*d_fd_dgd/d_fd_dfdd);
        fdd += diff_fdd;
        gd  += diff_gd;
        //println!(" iter {:?} err {:?} state_final= {:#?}", iterations, error, state_dgd);
        //println!(" fdd {:#?} change ({:#?}) gd {:#?} change ({:#?})", fdd, diff_fdd, gd, diff_gd);
        //println!(" fd_err {:#?} g_err {:#?}", fd_err, g_err);

        iterations += 1;
        if iterations>100 { panic!("Too many iterations of newton solve"); }
    }

    println!("Got fdd {:#?} gd {:#?}", fdd, gd);
    let fdd_final = Complex64::new(fdd, 0.0);
    let gd_final = Complex64::new(gd, 0.0);
    let state_final = integrate_through_bl(fdd_final, gd_final, &pm);

    println!("Final state {:#?}", state_final);

    //println!("Error in final values of fd= {:?} g= {:?}", state_final.fd.re-1.0,
    //                                                      state_final.g.re-1.0);
    //println!("d_g_dfdd:  {:#?}", d_g_dfdd);
    //println!("d_fd_dfdd: {:#?}", d_fd_dfdd);
    //println!("d_g_dgd:   {:#?}", d_g_dgd);
    //println!("d_fd_dgd:  {:#?}", d_fd_dgd);

    println!("Done.");
}
