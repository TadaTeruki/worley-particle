use super::{
    lerp::{InterpolationMethod, ParticleMapAttributeLerp},
    ParticleMap, ParticleMapAttribute,
};

pub trait ParticleMapAttributeGrad:
    ParticleMapAttribute + ParticleMapAttributeLerp + Default + PartialOrd
{
    fn grad(&self, other: &Self, delta: f64) -> Self;
}

impl ParticleMapAttributeGrad for f64 {
    fn grad(&self, other: &Self, delta: f64) -> Self {
        (other - self) / delta
    }
}

impl ParticleMapAttributeGrad for () {
    fn grad(&self, _: &Self, _: f64) -> Self {
        ()
    }
}

pub struct GradStrategy {
    pub delta: f64,
    pub iteration: usize,
    pub sample_num: usize,
}

pub struct Gradient<T: ParticleMapAttributeGrad> {
    pub angle: f64,
    pub value: T,
}

impl<T: ParticleMapAttributeGrad> ParticleMap<T> {
    pub fn as_unit_vector(&self, angle: f64) -> (f64, f64) {
        (angle.cos(), angle.sin())
    }
}

impl<T: ParticleMapAttributeGrad> ParticleMap<T> {
    pub fn get_forward_gradient(
        &self,
        x: f64,
        y: f64,
        strategy: &GradStrategy,
        method: &InterpolationMethod,
    ) -> Option<Gradient<T>> {
        let mut final_angle = 0.0;
        let mut final_value = T::default();
        let mut range = (0., std::f64::consts::PI * 2.);
        let center_value = self.get_interpolated(x, y, method)?;

        for _ in 0..strategy.iteration {
            let mut min_value = None;
            let mut min_angle = 0.0;
            let step = (range.1 - range.0) / (strategy.sample_num - 1) as f64;

            for i in 0..strategy.sample_num {
                let angle = range.0 + step * (i as f64);
                let dx = angle.cos() * strategy.delta;
                let dy = angle.sin() * strategy.delta;

                if let Some(value) = self.get_interpolated(x + dx, y + dy, method) {
                    if min_value.is_none() || &value < min_value.as_ref().unwrap() {
                        min_value = Some(value);
                        min_angle = angle;
                    }
                }
            }

            if min_value.is_none() {
                return None;
            }

            range = (min_angle - step * 0.5, min_angle + step * 0.5);
            final_angle = min_angle;
            final_value = min_value.unwrap();
        }

        let gradient_value = center_value.grad(&final_value, strategy.delta);

        Some(Gradient {
            angle: final_angle,
            value: gradient_value,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::lerp::InterpolationMethod;
    use super::*;
    use crate::{map::lerp::IDWStrategy, Particle, ParticleParameters};
    use std::collections::HashMap;

    #[test]
    fn test_get_gradient_simple_slope() {
        let params = ParticleParameters::new(0.0, 0.0, 1.0, 0).unwrap();
        let mut particles = HashMap::new();

        for y in -1..=1 {
            for x in -1..=1 {
                let particle = Particle::new(x, y, params);
                let value = x as f64;
                particles.insert(particle, value);
            }
        }

        let map = ParticleMap::new(params, particles);

        let strategy = GradStrategy {
            delta: 0.5,
            iteration: 3,
            sample_num: 8,
        };

        let gradient = map
            .get_forward_gradient(
                0.0,
                0.0,
                &strategy,
                &InterpolationMethod::IDW(IDWStrategy::default_from_params(&params)),
            )
            .unwrap();

        let expected_angle = std::f64::consts::PI;

        let gradient_unit = map.as_unit_vector(gradient.angle);
        let expected_unit = map.as_unit_vector(expected_angle);
        let unit_diff =
            (gradient_unit.0 - expected_unit.0).hypot(gradient_unit.1 - expected_unit.1);
        assert!(
            unit_diff < 0.1,
            "Expected angle near {:?}, got {:?}",
            expected_unit,
            gradient_unit
        );

        let expected_value = -1.0;

        assert!(
            (gradient.value - expected_value).abs() < 0.1,
            "Expected gradient value near {}, got {}",
            expected_value,
            gradient.value
        );
    }

    #[test]
    fn test_get_gradient_no_interpolation() {
        let params = ParticleParameters::new(0.0, 0.0, 1.0, 0).unwrap();
        let map: ParticleMap<f64> = ParticleMap::new(params, HashMap::new()); // Empty map

        let strategy = GradStrategy {
            delta: 0.1,
            iteration: 2,
            sample_num: 8,
        };

        let gradient = map.get_forward_gradient(0.0, 0.0, &strategy, &InterpolationMethod::Nearest);
        assert!(gradient.is_none(), "Expected None for empty map");
    }
}
