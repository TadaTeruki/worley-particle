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

impl<T: ParticleMapAttribute> ParticleMap<T> {
    pub fn new(params: ParticleParameters, particles: HashMap<Particle, T>) -> Self {
        Self { params, particles }
    }

    pub fn iter(&self) -> impl Iterator<Item = (&Particle, &T)> {
        self.particles.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&Particle, &mut T)> {
        self.particles.iter_mut()
    }

    pub fn insert(&mut self, particle: Particle, value: T) {
        self.particles.insert(particle, value);
    }

    pub fn remove(&mut self, particle: &Particle) -> Option<T> {
        self.particles.remove(particle)
    }

    pub fn get_mut(&mut self, particle: &Particle) -> Option<&mut T> {
        self.particles.get_mut(particle)
    }

    pub fn get_or_insert_with(&mut self, particle: Particle, f: impl FnOnce() -> T) -> &mut T {
        self.particles.entry(particle).or_insert_with(f)
    }

    pub fn get_or_insert(&mut self, particle: Particle, value: T) -> &mut T {
        self.particles.entry(particle).or_insert(value)
    }

    pub fn get_or_insert_default(&mut self, particle: Particle) -> &mut T
    where
        T: Default,
    {
        self.particles.entry(particle).or_insert_with(T::default)
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

        let second = second.iter().filter_map(|(particle, value)| {
            if self.particles.contains_key(particle) {
                None
            } else {
                Some((*particle, value.clone()))
            }
        });

        let map = self
            .particles
            .clone()
            .into_iter()
            .chain(second)
            .collect::<HashMap<_, _>>();

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
                        self.particles
                            .get(&particle)
                            .map(|value| (particle, value.clone()))
                    })
                    .collect::<HashMap<_, _>>();

                Self::new(self.params.clone(), particles)
            })
            .collect();

        maps
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
        let mut map1 = ParticleMap::new(params, HashMap::new());
        map1.insert(p1, "A");
        map1.insert(p2, "B");

        let p3 = Particle::new(0, 1, params);
        let p4 = Particle::new(1, 1, params);
        let mut map2 = ParticleMap::new(params, HashMap::new());
        map2.insert(p2, "X");
        map2.insert(p3, "C");
        map2.insert(p4, "D");

        let chained = map1.chain(&map2);
        assert!(chained.is_some());
        let chained = chained.unwrap();

        assert_eq!(chained.particles.len(), 4);
        assert_eq!(chained.get(&p1), Some(&"A"));
        assert_eq!(chained.get(&p2), Some(&"B"));
        assert_eq!(chained.get(&p3), Some(&"C"));
        assert_eq!(chained.get(&p4), Some(&"D"));

        let different_params = ParticleParameters::new(0.3, 0.3, 1.0, 0).unwrap();
        let mut map3 = ParticleMap::new(different_params, HashMap::new());
        map3.insert(Particle::new(0, 0, different_params), "E");

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

        let mut map = ParticleMap::new(params, HashMap::new());
        map.insert(p1, 1);
        map.insert(p2, 2);
        map.insert(p3, 3);
        map.insert(p4, 4);
        map.insert(p5, 5);

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
