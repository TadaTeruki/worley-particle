use crate::{Particle, ParticleParameters};

pub mod idw;
pub mod nni;

pub trait InterpolationStrategy {
    fn calculate_weights(
        &self,
        x: f64,
        y: f64,
        params: ParticleParameters,
    ) -> Option<Vec<(Particle, f64)>>;
}
