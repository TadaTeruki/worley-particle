use crate::{map::ParticleMap, map::ParticleMapAttribute, Particle, ParticleParameters};

/// The strategy for the Inverse Distance Weighting (IDW) interpolation.
#[derive(Debug, Clone, PartialEq, Copy)]
pub struct IDWStrategy {
    pub sample_max_distance: f64,
    pub weight_power: f64,
    pub smooth_power: Option<f64>,
    pub tolerance: f64,
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
}

enum IDWWeight {
    Inside(f64),
    Equal,
    Outside,
}

/// Calculate the IDW weight for a single point.
fn calculate_idw_weight(
    x: f64,
    y: f64,
    site: (f64, f64),
    sample_max_distance: f64,
    tolerance: f64,
    weight_power: f64,
    smooth_power: Option<f64>,
) -> IDWWeight {
    let distance = ((x - site.0).powi(2) + (y - site.1).powi(2)).sqrt();
    if distance > sample_max_distance {
        return IDWWeight::Outside;
    }
    if distance < tolerance {
        return IDWWeight::Equal;
    }
    let weight = if let Some(smooth_power) = smooth_power {
        (1.0 - (distance / sample_max_distance).powf(smooth_power)) / distance.powf(weight_power)
    } else {
        distance.powf(-weight_power)
    };

    IDWWeight::Inside(weight)
}

/// Calculate the IDW weights of the particles around the given site.
pub fn calculate_idw_weights<T: ParticleMapAttribute>(
    particle_map: &ParticleMap<T>,
    x: f64,
    y: f64,
    strategy: &IDWStrategy,
) -> Option<Vec<(Particle, f64)>> {
    let inside =
        Particle::from_inside_radius(x, y, *particle_map.params(), strategy.sample_max_distance);

    let particles_iter = inside
        .iter()
        .filter_map(|particle| particle_map.get(particle).map(|value| (*particle, value)));

    let mut weights = vec![];

    for (particle, _) in particles_iter {
        let weight = calculate_idw_weight(
            x,
            y,
            particle.site(),
            strategy.sample_max_distance,
            strategy.tolerance,
            strategy.weight_power,
            strategy.smooth_power,
        );

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
