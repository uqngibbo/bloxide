/*
    Code for computing viscosity with the Sutherland expression 

    @author: Nick Gibbons
*/

use num_complex::ComplexFloat;
use crate::mixedarithmetic::{Cplx,Mxd};

const MU_REF: f64 = 1.716e-05;
const T_REF: f64 = 273.0;
const S: f64 = 111.0;

// T here is a generic type, not the temperature!
pub fn sutherland_mu<T: ComplexFloat>(TEMP: T) -> T where T: Cplx<T>, f64: Mxd<T> {
    return MU_REF*ComplexFloat::sqrt(TEMP/T_REF)*(TEMP/T_REF)*(T_REF + S)/(TEMP + S);
}

pub fn sutherland_mu_derivative<T: ComplexFloat>(TEMP: T) -> T where T: Cplx<T>, f64: Mxd<T> {
    return MU_REF*(T_REF+S)*ComplexFloat::sqrt(TEMP/T_REF)*(3.0*S+TEMP)/(2.0*T_REF*(S+TEMP)*(S+TEMP));
}

