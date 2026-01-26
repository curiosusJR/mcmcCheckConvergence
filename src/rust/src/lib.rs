#![allow(nonstandard_style)]

pub mod core;
pub mod tree;

#[cfg(feature = "r")]
mod r_api;
#[cfg(feature = "r")]
pub use r_api::*;
