/*
    Let's use structs with named fields instead of array. I can't 
    get the complex numbers to work properly with those.

    @author: Nick Gibbons
*/

use num_complex::Complex64;
use std::ops::{Add,Sub,Div,Mul};

// TODO: Maybe we don't want copy here???
#[derive(PartialEq, Debug, Clone, Copy)]
pub struct State {
    pub f: Complex64,
    pub fd: Complex64,
    pub fdd: Complex64,
    pub g: Complex64,
    pub gd: Complex64,
    pub y: Complex64,
}

impl Add for State{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            f:   self.f + rhs.f,
            fd:  self.fd + rhs.fd,
            fdd: self.fdd + rhs.fdd,
            g:   self.g + rhs.g,
            gd:  self.gd + rhs.gd,
            y:   self.y + rhs.y,
        }
    }
}

impl Add<Complex64> for State{
    type Output = Self;
    fn add(self, rhs: Complex64) -> Self {
        Self {
            f:   self.f + rhs,
            fd:  self.fd + rhs,
            fdd: self.fdd + rhs,
            g:   self.g + rhs,
            gd:  self.gd + rhs,
            y:    self.y + rhs,
        }
    }
}

impl Add<f64> for State{
    type Output = Self;
    fn add(self, rhs: f64) -> Self {
        Self {
            f:   self.f + rhs,
            fd:  self.fd + rhs,
            fdd: self.fdd + rhs,
            g:   self.g + rhs,
            gd:  self.gd + rhs,
            y:    self.y + rhs,
        }
    }
}

impl Mul<State> for State{
    type Output = Self;
    fn mul(self, rhs: State) -> Self {
        Self {
            f:   self.f*rhs.f,
            fd:  self.fd*rhs.fd,
            fdd: self.fdd*rhs.fdd,
            g:   self.g*rhs.g,
            gd:  self.gd*rhs.gd,
            y:   self.y*rhs.y,
        }
    }
}

impl Mul<Complex64> for State{
    type Output = Self;
    fn mul(self, rhs: Complex64) -> Self {
        Self {
            f:   self.f*rhs,
            fd:  self.fd*rhs,
            fdd: self.fdd*rhs,
            g:   self.g*rhs,
            gd:  self.gd*rhs,
            y:   self.y*rhs,
        }
    }
}

impl Mul<State> for Complex64 {
    type Output = State;
    fn mul(self, rhs: State) -> State {
        State {
            f:   self*rhs.f,
            fd:  self*rhs.fd,
            fdd: self*rhs.fdd,
            g:   self*rhs.g,
            gd:  self*rhs.gd,
            y:   self*rhs.y,
        }
    }
}

impl Mul<f64> for State{
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            f:   self.f*rhs,
            fd:  self.fd*rhs,
            fdd: self.fdd*rhs,
            g:   self.g*rhs,
            gd:  self.gd*rhs,
            y:   self.y*rhs,
        }
    }
}

impl Mul<State> for f64 {
    type Output = State;
    fn mul(self, rhs: State) -> State {
        State {
            f:   self*rhs.f,
            fd:  self*rhs.fd,
            fdd: self*rhs.fdd,
            g:   self*rhs.g,
            gd:  self*rhs.gd,
            y:   self*rhs.y,
        }
    }
}

impl Sub for State{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            f:   self.f  - rhs.f,
            fd:  self.fd - rhs.fd,
            fdd: self.fdd- rhs.fdd,
            g:   self.g  - rhs.g,
            gd:  self.gd - rhs.gd,
            y:   self.y  - rhs.y,
        }
    }
}

impl Div for State{
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        Self {
            f:   self.f/rhs.f,
            fd:  self.fd/rhs.fd,
            fdd: self.fdd/rhs.fdd,
            g:   self.g/rhs.g,
            gd:  self.gd/rhs.gd,
            y:   self.y/rhs.y,
        }
    }
}

impl Div<State> for Complex64 {
    type Output = State;
    fn div(self, rhs: State) -> State {
        State {
            f:   self/rhs.f,
            fd:  self/rhs.fd,
            fdd: self/rhs.fdd,
            g:   self/rhs.g,
            gd:  self/rhs.gd,
            y:   self/rhs.y,
        }
    }
}

impl Div<Complex64> for State{
    type Output = Self;
    fn div(self, rhs: Complex64) -> Self {
        Self {
            f:   self.f/rhs,
            fd:  self.fd/rhs,
            fdd: self.fdd/rhs,
            g:   self.g/rhs,
            gd:  self.gd/rhs,
            y:   self.y/rhs,
        }
    }
}

impl Div<State> for f64 {
    type Output = State;
    fn div(self, rhs: State) -> State {
        State {
            f:   self/rhs.f,
            fd:  self/rhs.fd,
            fdd: self/rhs.fdd,
            g:   self/rhs.g,
            gd:  self/rhs.gd,
            y:   self/rhs.y,
        }
    }
}

impl Div<f64> for State{
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self {
            f:   self.f/rhs,
            fd:  self.fd/rhs,
            fdd: self.fdd/rhs,
            g:   self.g/rhs,
            gd:  self.gd/rhs,
            y:   self.y/rhs,
        }
    }
}

impl State{
    pub fn abs(self) -> Self {
        Self {
            f:   Complex64::new(f64::abs(self.f.re),   self.f.im),
            fd:  Complex64::new(f64::abs(self.fd.re),  self.fd.im),
            fdd: Complex64::new(f64::abs(self.fdd.re), self.fdd.im),
            g:   Complex64::new(f64::abs(self.g.re),   self.g.im),
            gd:  Complex64::new(f64::abs(self.gd.re),  self.gd.im),
            y:   Complex64::new(f64::abs(self.y.re),   self.y.im),
        }
    }
}

impl State{
    pub fn wall_state(fdd: f64, gd: f64, h_wall: f64, h_e: f64 ) -> Self {
        Self {
            f:   Complex64::new(0.0, 0.0),
            fd:  Complex64::new(0.0, 0.0),
            fdd: Complex64::new(fdd, 0.0),
            g:   Complex64::new(h_wall/h_e, 0.0),
            gd:  Complex64::new(gd, 0.0),
            y:   Complex64::new(0.0, 0.0),
        }
    }
}

