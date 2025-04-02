use crate::{Particle, ParticleParameters};

/// The strategy for the Inverse Distance Weighting (IDW) interpolation.
#[derive(Debug, Clone, PartialEq, Copy)]
pub struct IDWStrategy {
    pub sample_max_distance: f64,
    pub weight_power: f64,
    pub smooth_power: Option<f64>,
    pub tolerance: f64,
}

enum IDWWeight {
    Inside(f64),
    Equal,
    Outside,
}

impl IDWStrategy {
    pub fn default_from_params(params: &ParticleParameters) -> Self {
        Self {
            sample_max_distance: params.scale * 1.415,
            tolerance: params.scale * 1e-6,
            weight_power: 1.5,
            smooth_power: Some(1.0),
        }
    }

    /// Calculate the IDW weight for a single point.
    fn calculate_idw_weight(&self, x: f64, y: f64, site: (f64, f64)) -> IDWWeight {
        let distance = ((x - site.0).powi(2) + (y - site.1).powi(2)).sqrt();
        if distance > self.sample_max_distance {
            return IDWWeight::Outside;
        }
        if distance < self.tolerance {
            return IDWWeight::Equal;
        }
        let weight = if let Some(smooth_power) = self.smooth_power {
            (1.0 - (distance / self.sample_max_distance).powf(smooth_power))
                / distance.powf(self.weight_power)
        } else {
            distance.powf(-self.weight_power)
        };

        IDWWeight::Inside(weight)
    }

    /// Calculate the IDW weights of the particles around the given site.
    pub fn calculate_idw_weights(
        &self,
        x: f64,
        y: f64,
        params: ParticleParameters,
    ) -> Option<Vec<(Particle, f64)>> {
        let inside = Particle::from_inside_radius(x, y, params, self.sample_max_distance);

        let mut weights = vec![];

        for particle in inside {
            let weight = self.calculate_idw_weight(x, y, particle.site());

            match weight {
                IDWWeight::Inside(w) => {
                    weights.push((particle, w));
                }
                IDWWeight::Equal => {
                    return Some(vec![(particle, 1.0)]);
                }
                IDWWeight::Outside => {}
            }
        }

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
