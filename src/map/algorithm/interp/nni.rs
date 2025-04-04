use naturalneighbor::{Interpolator, Point};

use crate::{
    map::{ParticleMap, ParticleMapAttribute},
    Particle, ParticleParameters,
};

use super::InterpolationStrategy;

/// The strategy for the Natural Neighbor Interpolation (NNI).
#[derive(Clone)]
pub struct NNIStrategy {
    kind: NNIStrategyKind,
}

#[derive(Clone)]
enum NNIStrategyKind {
    Prebuilt {
        interpolator: Interpolator,
        particles: Vec<Particle>,
    },
    Immediate,
}

impl NNIStrategy {
    pub fn new_immediate() -> Self {
        Self {
            kind: NNIStrategyKind::Immediate,
        }
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

        Self {
            kind: NNIStrategyKind::Prebuilt {
                interpolator: Interpolator::new(&points),
                particles,
            },
        }
    }
}

impl InterpolationStrategy for NNIStrategy {
    fn calculate_weights(
        &self,
        x: f64,
        y: f64,
        params: ParticleParameters,
    ) -> Option<Vec<(Particle, f64)>> {
        let particles = match &self.kind {
            NNIStrategyKind::Prebuilt { particles, .. } => particles,
            NNIStrategyKind::Immediate => {
                &Particle::from_inside_radius(x, y, params, params.scale * 2.0)
            }
        };

        let interpolator = match &self.kind {
            NNIStrategyKind::Prebuilt { interpolator, .. } => interpolator,
            NNIStrategyKind::Immediate => {
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
