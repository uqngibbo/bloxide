/*
    Storage area for horrible mixed arithmetic traits that placate the rust compiler
    into letting me do mixed complex/real arithmetic

    @author: Nick Gibbons
*/

use std::ops::{Add,Mul,Div,Sub};

// Super traits to allow mixed complex/real arithmetic
pub trait Cplx<T>: Mul<f64, Output = T> + Add<f64, Output = T>   + Div<f64, Output = T>   + Sub<f64, Output = T> {}
impl<T> Cplx<T> for T where T: Mul<f64, Output = T>   + Add<f64, Output = T>   + Div<f64, Output = T>   + Sub<f64, Output = T> {}


pub trait Mxd<T>:    Mul<T,Output = T>   + Add<T,   Output = T>   + Div<T,   Output = T>   + Sub<T,   Output = T>
               + Mul<f64, Output = f64> + Add<f64, Output = f64> + Div<f64, Output = f64> + Sub<f64, Output = f64> {}
impl<T> Mxd<T> for f64 where f64: Mul<T, Output = T>   + Add<T, Output = T>   + Div<T, Output = T>   + Sub<T, Output = T>
                                + Mul<f64, Output = f64> + Add<f64, Output = f64> + Div<f64, Output = f64> + Sub<f64, Output = f64> {}
