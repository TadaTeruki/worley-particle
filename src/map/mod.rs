use algorithm::disjoint_set::DisjointSet;

use crate::{Particle, ParticleParameters};
use std::{collections::HashMap, fmt::Debug};

mod algorithm;
pub mod grad;
pub mod lerp;
pub mod network;
pub mod rw;

pub trait ParticleMapAttribute: Debug + Clone + PartialEq {}
impl<T: Debug + Clone + PartialEq> ParticleMapAttribute for T {}

#[derive(Debug, Clone)]
pub struct ParticleMap<T: ParticleMapAttribute> {
    params: ParticleParameters,
    particles: HashMap<Particle, T>,
}

impl<T: ParticleMapAttribute> FromIterator<(Particle, T)> for ParticleMap<T> {
    fn from_iter<I: IntoIterator<Item = (Particle, T)>>(iter: I) -> Self {
        let particles = iter.into_iter().collect::<HashMap<_, _>>();
        let params = particles
            .keys()
            .next()
            .map_or_else(|| ParticleParameters::default(), |particle| particle.params);
        Self { params, particles }
    }
}

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

impl<T: ParticleMapAttribute> ParticleMap<T> {
    /// Create a new particle map with the given parameters and particles.
    ///
    /// Note that the particles with different parameters will be ignored.
    pub fn new(params: ParticleParameters, particles: HashMap<Particle, T>) -> Self {
        Self {
            params,
            particles: particles
                .into_iter()
                .filter(|(particle, _)| particle.params == params)
                .collect(),
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (&Particle, &T)> {
        self.particles.iter()
    }

    /// Returns an iterator that clones the keys and values.
    pub fn cloned_iter(&self) -> impl Iterator<Item = (Particle, T)> + '_ {
        self.particles.iter().map(|(k, v)| (*k, v.clone()))
    }

    /// Get the sites of the particles in the map.
    pub fn sites(&self) -> Vec<(f64, f64)> {
        self.particles
            .keys()
            .map(|particle| particle.site())
            .collect()
    }

    /// Get the parameters of the particle map.
    pub fn params(&self) -> &ParticleParameters {
        &self.params
    }

    /// Get the value of the particle at the given site.
    pub fn get(&self, particle: &Particle) -> Option<&T> {
        self.particles.get(particle)
    }

    /// Get corner coordinates of the range of the registered particles.
    pub fn corners(&self) -> ((f64, f64), (f64, f64)) {
        let particle_sites = self
            .particles
            .keys()
            .map(|particle| particle.site())
            .collect::<Vec<_>>();

        let min_x = particle_sites
            .iter()
            .map(|(x, _)| *x)
            .fold(f64::MAX, |acc, x| acc.min(x));

        let max_x = particle_sites
            .iter()
            .map(|(x, _)| *x)
            .fold(f64::MIN, |acc, x| acc.max(x));

        let min_y = particle_sites
            .iter()
            .map(|(_, y)| *y)
            .fold(f64::MAX, |acc, y| acc.min(y));

        let max_y = particle_sites
            .iter()
            .map(|(_, y)| *y)
            .fold(f64::MIN, |acc, y| acc.max(y));

        ((min_x, min_y), (max_x, max_y))
    }

    /// Chain the map with another map if the parameters are same.
    /// The particles in the second map which are already in the first map will be ignored.
    pub fn chain(&self, second: &Self) -> Option<Self> {
        if self.params != second.params {
            return None;
        }

        let second = second
            .cloned_iter()
            .filter(|(particle, _)| !self.particles.contains_key(particle));

        let map = self.cloned_iter().chain(second).collect::<HashMap<_, _>>();

        Some(Self::new(self.params.clone(), map))
    }

    /// Divide the map into clusters of connected particles.
    pub fn divide_into_clusters(&self) -> Vec<Self> {
        let keys = self.particles.keys().cloned().collect::<Vec<_>>();
        let mut disjoint_set = DisjointSet::new(&keys.clone());

        for &particle in &keys {
            for neighbor in particle.calculate_voronoi().neighbors {
                if self.particles.contains_key(&neighbor) {
                    disjoint_set.union(particle, neighbor);
                }
            }
        }

        let sets = disjoint_set.get_all_sets();

        let maps = sets
            .into_iter()
            .map(|set| {
                let particles = set
                    .into_iter()
                    .filter_map(|particle| {
                        if let Some(value) = self.particles.get(&particle) {
                            Some((particle, value.clone()))
                        } else {
                            None
                        }
                    })
                    .collect::<HashMap<_, _>>();

                Self::new(self.params.clone(), particles)
            })
            .collect();

        maps
    }

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
            (1.0 - (distance / sample_max_distance).powf(smooth_power))
                / distance.powf(weight_power)
        } else {
            distance.powf(-weight_power)
        };

        IDWWeight::Inside(weight)
    }

    /// Calculate the IDW weights of the particles around the given site.
    pub fn calculate_idw_weights(
        &self,
        x: f64,
        y: f64,
        strategy: &IDWStrategy,
    ) -> Option<Vec<(Particle, f64)>> {
        let inside = Particle::from_inside_square(x, y, self.params, strategy.sample_max_distance);
        let particles_iter = inside
            .iter()
            .filter_map(|particle| self.particles.get(particle).map(|value| (*particle, value)));

        let mut weights = vec![];

        for (particle, _) in particles_iter {
            let weight = Self::calculate_idw_weight(
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

    #[allow(dead_code)]
    fn is_same(&self, other: &Self) -> bool {
        if self.particles.len() != other.particles.len() {
            return false;
        }

        for (particle, value) in self.particles.iter() {
            if let Some(other_value) = other.particles.get(particle) {
                if value != other_value {
                    return false;
                }
            } else {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Particle;

    #[test]
    fn test_chain() {
        let params = ParticleParameters::default();

        let p1 = Particle::new(0, 0, params);
        let p2 = Particle::new(1, 0, params);
        let hash_map = vec![(p1, "A"), (p2, "B")].into_iter().collect();
        let map1 = ParticleMap::new(params, hash_map);

        let p3 = Particle::new(0, 1, params);
        let p4 = Particle::new(1, 1, params);
        let hash_map = vec![(p2, "X"), (p3, "C"), (p4, "D")].into_iter().collect();
        let map2 = ParticleMap::new(params, hash_map);

        let chained = map1.chain(&map2);
        assert!(chained.is_some());
        let chained = chained.unwrap();

        assert_eq!(chained.particles.len(), 4);
        assert_eq!(chained.get(&p1), Some(&"A"));
        assert_eq!(chained.get(&p2), Some(&"B"));
        assert_eq!(chained.get(&p3), Some(&"C"));
        assert_eq!(chained.get(&p4), Some(&"D"));

        let different_params = ParticleParameters::new(0.3, 0.3, 1.0, 0).unwrap();
        let hash_map = vec![(Particle::new(0, 0, different_params), "E")]
            .into_iter()
            .collect();

        let map3 = ParticleMap::new(different_params, hash_map);

        let chained = map1.chain(&map3);
        assert!(chained.is_none());
    }

    #[test]
    fn test_divide_into_clusters() {
        let params = ParticleParameters::default();

        let p1 = Particle::new(0, 0, params);
        let p2 = Particle::new(1, 0, params);
        let p3 = Particle::new(0, 1, params);
        let p4 = Particle::new(3, 3, params);
        let p5 = Particle::new(4, 3, params);

        let hash_map = vec![(p1, 1), (p2, 2), (p3, 3), (p4, 4), (p5, 5)]
            .into_iter()
            .collect();

        let map = ParticleMap::new(params, hash_map);

        let clusters = map.divide_into_clusters();

        assert_eq!(clusters.len(), 2);

        let cluster_sizes: Vec<usize> = clusters
            .iter()
            .map(|cluster| cluster.particles.len())
            .collect();
        assert!(cluster_sizes.contains(&3));
        assert!(cluster_sizes.contains(&2));

        for cluster in &clusters {
            if cluster.particles.len() == 3 {
                assert!(cluster.get(&p1).is_some());
                assert!(cluster.get(&p2).is_some());
                assert!(cluster.get(&p3).is_some());
            } else if cluster.particles.len() == 2 {
                assert!(cluster.get(&p4).is_some());
                assert!(cluster.get(&p5).is_some());
            }
        }
    }
}
