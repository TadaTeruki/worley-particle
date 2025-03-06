use crate::{Particle, ParticleParameters};
use std::{collections::HashMap, fmt::Debug};

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

    pub fn sites(&self) -> Vec<(f64, f64)> {
        self.particles
            .keys()
            .map(|particle| particle.site())
            .collect()
    }

    pub fn params(&self) -> &ParticleParameters {
        &self.params
    }

    pub fn get(&self, particle: &Particle) -> Option<&T> {
        self.particles.get(particle)
    }

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
