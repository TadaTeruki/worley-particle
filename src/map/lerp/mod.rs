use std::collections::HashMap;

use crate::{Particle, ParticleParameters};

use super::{
    algorithm::interp::{idw::IDWStrategy, nni::NNIStrategy, InterpolationStrategy},
    ParticleMap, ParticleMapAttribute,
};

pub mod vertorization;

pub trait ParticleMapAttributeLerp: ParticleMapAttribute + Default {
    fn lerp(&self, other: &Self, t: f64) -> Self;
}

impl ParticleMapAttributeLerp for f64 {
    fn lerp(&self, other: &Self, t: f64) -> Self {
        self + (other - self) * t
    }
}

impl ParticleMapAttributeLerp for () {
    fn lerp(&self, _: &Self, _: f64) -> Self {}
}

pub enum InterpolationMethod {
    Nearest,
    IDW(IDWStrategy),
    IDWSeparated(IDWStrategy),
    NNI(NNIStrategy),
}

impl<T: ParticleMapAttributeLerp> ParticleMap<T> {
    pub fn get_interpolated(&self, x: f64, y: f64, method: &InterpolationMethod) -> Option<T> {
        match method {
            InterpolationMethod::Nearest => {
                let particle = Particle::from(x, y, self.params);
                self.particles.get(&particle).cloned()
            }
            InterpolationMethod::IDW(strategy) => {
                let weights = strategy.calculate_weights(x, y, self.params)?;
                let mut total_value: Option<T> = None;
                let mut tmp_weight = 0.0;

                for (particle, weight) in weights {
                    let value = if let Some(value) = self.particles.get(&particle) {
                        value
                    } else {
                        continue;
                    };
                    tmp_weight += weight;

                    if let Some(total_value) = total_value.as_mut() {
                        *total_value = total_value.lerp(value, weight / tmp_weight);
                    } else {
                        total_value = Some(value.clone());
                    }
                }

                total_value
            }
            InterpolationMethod::NNI(strategy) => {
                let weights = strategy.calculate_weights(x, y, self.params)?;
                let mut total_value: Option<T> = None;
                let mut tmp_weight = 0.0;

                for (particle, weight) in weights {
                    let value = if let Some(value) = self.particles.get(&particle) {
                        value
                    } else {
                        continue;
                    };
                    tmp_weight += weight;

                    if let Some(total_value) = total_value.as_mut() {
                        *total_value = total_value.lerp(value, weight / tmp_weight);
                    } else {
                        total_value = Some(value.clone());
                    }
                }

                total_value
            }
            InterpolationMethod::IDWSeparated(strategy) => {
                let weights = strategy.calculate_weights(x, y, self.params)?;
                let self_particle = Particle::from(x, y, self.params);

                let weight_sum = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
                let mut total_value: Option<T> = None;
                let mut tmp_weight = 0.0;
                for (particle, weight) in weights {
                    let (value, weight) = if particle == self_particle {
                        (
                            self.particles.get(&particle).unwrap(),
                            (weight / weight_sum - 0.5).max(0.0),
                        )
                    } else {
                        (&T::default(), weight / weight_sum)
                    };
                    tmp_weight += weight;
                    if tmp_weight <= 0.0 {
                        continue;
                    }
                    if let Some(total_value) = total_value.as_mut() {
                        *total_value = total_value.lerp(value, weight / tmp_weight);
                    } else {
                        total_value = Some(value.clone());
                    }
                }

                total_value
            }
        }
    }

    /// Chain the map with another map.
    /// If the parameters are different, the second map will be interpolated to the first map.
    /// The particles in the second map which are already in the first map will be ignored.
    pub fn force_chain_with_interpolation(
        &self,
        second: &Self,
        interp_method: InterpolationMethod,
    ) -> Self {
        let second = if self.params == second.params {
            second
        } else {
            &second
                .map_with_another_params_iter(interp_method, self.params)
                .collect::<ParticleMap<_>>()
        };

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

        Self::new(self.params, map)
    }

    /// Chain the map with another map if the parameters are same.
    /// If the parameters are different, the second map will be interpolated with nearest method to the first map.
    /// The particles in the second map which are already in the first map will be ignored.
    pub fn force_chain(&self, second: &Self) -> Self {
        self.force_chain_with_interpolation(second, InterpolationMethod::Nearest)
    }

    pub fn map_with_another_params_iter(
        &self,
        interp_method: InterpolationMethod,
        another_params: ParticleParameters,
    ) -> impl Iterator<Item = (Particle, T)> + '_ {
        self.particles
            .iter()
            .flat_map(move |(&range_particle, _)| {
                Particle::from_inside_particle(another_params, range_particle)
            })
            .filter_map(move |particle| {
                let (x, y) = particle.site();
                let value = self.get_interpolated(x, y, &interp_method)?;
                Some((particle, value))
            })
    }

    pub fn rasterise(
        &self,
        img_width: usize,
        img_height: usize,
        corners: ((f64, f64), (f64, f64)),
        interp_method: InterpolationMethod,
    ) -> Vec<Vec<Option<T>>> {
        let ((mut min_x, mut min_y), (mut max_x, mut max_y)) = corners;
        if min_x > max_x {
            std::mem::swap(&mut min_x, &mut max_x);
        }
        if min_y > max_y {
            std::mem::swap(&mut min_y, &mut max_y);
        }
        let mut raster = vec![vec![None; img_width]; img_height];

        for (iy, item) in raster.iter_mut().enumerate().take(img_height) {
            for (ix, item) in item.iter_mut().enumerate().take(img_width) {
                let x = min_x + (max_x - min_x) * ix as f64 / img_width as f64;
                let y = min_y + (max_y - min_y) * iy as f64 / img_height as f64;
                *item = self.get_interpolated(x, y, &interp_method);
            }
        }

        raster
    }
}
