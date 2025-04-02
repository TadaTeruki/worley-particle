use naturalneighbor::{Interpolator, Point};

use crate::{
    map::{ParticleMap, ParticleMapAttribute},
    Particle, ParticleParameters,
};

/// The strategy for the Natural Neighbor Interpolation (NNI).
#[derive(Clone)]
pub enum NNIStrategy {
    Prebuilt {
        interpolator: Interpolator,
        particles: Vec<Particle>,
    },
    Immediate,
}

impl NNIStrategy {
    pub fn new_immediate() -> Self {
        Self::Immediate
    }

    pub fn new_prebuild<T: ParticleMapAttribute>(particle_map: &ParticleMap<T>) -> Self {
        let particles = particle_map.particles.keys().cloned().collect::<Vec<_>>();

        let points = particles
            .iter()
            .map(|particle| {
                let site = particle.site();
                Point {
                    x: site.0,
                    y: site.1,
                }
            })
            .collect::<Vec<_>>();

        Self::Prebuilt {
            interpolator: Interpolator::new(&points),
            particles,
        }
    }

    /// Calculate the NNI weights of the particles around the given site.
    pub fn calculate_nni_weights(
        &self,
        x: f64,
        y: f64,
        params: ParticleParameters,
    ) -> Option<Vec<(Particle, f64)>> {
        let particles = match self {
            Self::Prebuilt { particles, .. } => particles,
            Self::Immediate => &Particle::from_inside_radius(x, y, params, params.scale * 2.0),
        };

        let interpolator = match self {
            Self::Prebuilt { interpolator, .. } => interpolator,
            Self::Immediate => {
                let points = particles
                    .iter()
                    .map(|particle| {
                        let site = particle.site();
                        Point {
                            x: site.0,
                            y: site.1,
                        }
                    })
                    .collect::<Vec<_>>();

                &Interpolator::new(&points)
            }
        };

        let mut weights = interpolator
            .query_weights(Point { x, y })
            .ok()??
            .iter()
            .map(|(idx, weight)| (particles[*idx], *weight))
            .collect::<Vec<_>>();

        let sum = weights.iter().map(|(_, w)| w).sum::<f64>();

        if sum == 0.0 {
            return None;
        }

        for (_, w) in &mut weights {
            *w /= sum;
        }

        if weights.is_empty() {
            None
        } else {
            Some(weights)
        }
    }
}
