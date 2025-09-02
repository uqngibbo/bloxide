/*
    Let's use structs with named fields instead of array. I can't 
    get the complex numbers to work properly with those.

    @author: Nick Gibbons
*/

use std::ops::{Add,Sub,Div,Mul};

// TODO: Maybe we don't want copy here???
#[derive(PartialEq, Debug, Clone, Copy)]
pub struct State {
    pub f: f64,
    pub fd: f64,
    pub fdd: f64,
    pub g: f64,
    pub gd: f64,
    pub y: f64,
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
            f:   f64::abs(self.f),
            fd:  f64::abs(self.fd),
            fdd: f64::abs(self.fdd),
            g:   f64::abs(self.g),
            gd:  f64::abs(self.gd),
            y:   f64::abs(self.y),
        }
    }
}
