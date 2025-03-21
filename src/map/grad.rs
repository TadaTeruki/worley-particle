use std::ops::{Neg, Sub};

use super::{
    lerp::{InterpolationMethod, ParticleMapAttributeLerp},
    ParticleMap, ParticleMapAttribute,
};

pub trait ParticleMapAttributeGrad:
    ParticleMapAttribute
    + ParticleMapAttributeLerp
    + Default
    + PartialOrd
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Copy
{
    fn div_f64(&self, other: f64) -> Self;
    fn abs(&self) -> Self;
}

impl ParticleMapAttributeGrad for f64 {
    fn div_f64(&self, other: f64) -> Self {
        *self / other
    }

    fn abs(&self) -> Self {
        f64::abs(*self)
    }
}

#[derive(Debug, Clone, Copy)]
pub enum GradDifferenceType {
    Central,
    Forward,
}

#[derive(Debug, Clone, Copy)]
pub enum GradDirectionType {
    ToDown,
    ToUp,
    Steepest,
}

/// A strategy for gradient calculation.
///
/// The time complexity of tasks that involve gradient calculation is O(n * m), where n is the number of iterations and m is the number of samples.
///
/// # Fields
///
/// * `delta` - The step size used in the gradient calculation.
/// * `difference_type` - The type of difference method used (e.g., forward, backward, central).
/// * `direction_type` - The type of direction method used (e.g., coordinate, random).
/// * `iteration` - The number of iterations to perform.
/// * `sample_num` - The number of samples to use in the gradient calculation. Note that the interval between samples for the central difference method is PI / sample_num, and for the forward difference method it is 2 * PI / sample_num (larger).
pub struct GradStrategy {
    pub delta: f64,
    pub difference_type: GradDifferenceType,
    pub direction_type: GradDirectionType,
    pub iteration: usize,
    pub sample_num: usize,
}

impl Default for GradStrategy {
    fn default() -> Self {
        return Self {
            delta: 1e-6,
            direction_type: GradDirectionType::ToDown,
            difference_type: GradDifferenceType::Central,
            iteration: 4,
            sample_num: 8,
        };
    }
}

/// A struct representing a gradient with a specific angle and value.
///
/// # Type Parameters
///
/// * `T` - A type that implements the `ParticleMapAttributeGrad` trait.
///
/// # Fields
///
/// * `angle` - The angle of the gradient in radians.
/// * `value` - The value associated with the gradient, of type `T`.
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
    pub fn get_gradient(
        &self,
        x: f64,
        y: f64,
        strategy: &GradStrategy,
        method: &InterpolationMethod,
    ) -> Option<Gradient<T>> {
        match strategy.difference_type {
            GradDifferenceType::Forward => self.get_gradient_forward(x, y, strategy, method),
            GradDifferenceType::Central => self.get_gradient_central(x, y, strategy, method),
        }
    }

    fn get_gradient_forward(
        &self,
        x: f64,
        y: f64,
        strategy: &GradStrategy,
        method: &InterpolationMethod,
    ) -> Option<Gradient<T>> {
        let center_value = self.get_interpolated(x, y, method)?;
        let mut final_angle = 0.0;
        let mut final_value = T::default();
        let mut range = (0., std::f64::consts::PI * 2.);

        for _ in 0..strategy.iteration {
            let mut extreme_value = None;
            let mut extreme_angle = 0.0;
            let step = (range.1 - range.0) / (strategy.sample_num - 1) as f64;

            for i in 0..strategy.sample_num {
                let angle = range.0 + step * (i as f64);
                let dx = angle.cos() * strategy.delta;
                let dy = angle.sin() * strategy.delta;

                if let Some(value) = self.get_interpolated(x + dx, y + dy, method) {
                    match strategy.direction_type {
                        GradDirectionType::ToDown => {
                            if extreme_value.is_none() || &value < extreme_value.as_ref().unwrap() {
                                extreme_value = Some(value - center_value);
                                extreme_angle = angle;
                            }
                        }
                        GradDirectionType::ToUp => {
                            if extreme_value.is_none() || &value > extreme_value.as_ref().unwrap() {
                                extreme_value = Some(value - center_value);
                                extreme_angle = angle;
                            }
                        }
                        GradDirectionType::Steepest => {
                            let diff = value - center_value;
                            if extreme_value.is_none()
                                || diff.abs() > extreme_value.as_ref().unwrap().abs()
                            {
                                extreme_value = Some(diff);
                                extreme_angle = angle;
                            }
                        }
                    }
                }
            }

            if extreme_value.is_none() {
                return None;
            }

            range = (extreme_angle - step * 0.5, extreme_angle + step * 0.5);
            final_angle = extreme_angle;
            final_value = extreme_value.unwrap();
        }

        let gradient_value = final_value.div_f64(strategy.delta);

        Some(Gradient {
            angle: final_angle,
            value: gradient_value,
        })
    }

    fn get_gradient_central(
        &self,
        x: f64,
        y: f64,
        strategy: &GradStrategy,
        method: &InterpolationMethod,
    ) -> Option<Gradient<T>> {
        let mut final_angle = 0.0;
        let mut final_value = T::default();
        let mut range = (0., std::f64::consts::PI);

        for _ in 0..strategy.iteration {
            let mut extreme_value = None;
            let mut extreme_angle = 0.0;
            let step = (range.1 - range.0) / (strategy.sample_num - 1) as f64;

            for i in 0..strategy.sample_num {
                let angle = range.0 + step * (i as f64);
                let dx = angle.cos() * strategy.delta * 0.5;
                let dy = angle.sin() * strategy.delta * 0.5;

                // Get values on both sides for central difference
                let forward_value = self.get_interpolated(x + dx, y + dy, method);
                let backward_value = self.get_interpolated(x - dx, y - dy, method);

                if let (Some(forward), Some(backward)) = (forward_value, backward_value) {
                    // Calculate central difference
                    let diff = forward - backward;

                    match strategy.direction_type {
                        GradDirectionType::ToDown => {
                            if extreme_value.is_none() || &diff < extreme_value.as_ref().unwrap() {
                                extreme_value = Some(diff);
                                extreme_angle = angle;
                            }
                            if extreme_value.is_none() || &(-diff) < extreme_value.as_ref().unwrap()
                            {
                                extreme_value = Some(-diff);
                                extreme_angle = angle + std::f64::consts::PI;
                            }
                        }
                        GradDirectionType::ToUp => {
                            if extreme_value.is_none() || &diff > extreme_value.as_ref().unwrap() {
                                extreme_value = Some(diff);
                                extreme_angle = angle;
                            }
                            if extreme_value.is_none() || &(-diff) > extreme_value.as_ref().unwrap()
                            {
                                extreme_value = Some(-diff);
                                extreme_angle = angle + std::f64::consts::PI;
                            }
                        }
                        GradDirectionType::Steepest => {
                            if extreme_value.is_none()
                                || diff.abs() > extreme_value.as_ref().unwrap().abs()
                            {
                                extreme_value = Some(diff);
                                extreme_angle = angle;
                            }
                        }
                    }
                }
            }

            if extreme_value.is_none() {
                return None;
            }

            range = (extreme_angle - step * 0.5, extreme_angle + step * 0.5);
            final_angle = extreme_angle;
            final_value = extreme_value.unwrap();
        }

        Some(Gradient {
            angle: final_angle,
            value: final_value.div_f64(strategy.delta),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::lerp::InterpolationMethod;
    use super::*;
    use crate::{map::IDWStrategy, Particle, ParticleParameters};
    use std::collections::HashMap;

    /// Helper struct for gradient test cases
    struct GradTestCase {
        /// Strategy configuration
        strategy: GradStrategy,
        /// Name of the test case for error messages
        name: &'static str,
        /// Expected angle(s) in radians
        expected_angles: Vec<f64>,
        /// Threshold for gradient value comparison
        value_threshold: f64,
        /// Whether gradient value should be greater than threshold
        greater_than_threshold: bool,
        /// Whether to check absolute value
        check_abs: bool,
    }

    #[test]
    fn test_get_gradient_simple_slope() {
        let params = ParticleParameters::new(0.0, 0.0, 1.0, 0).unwrap();
        let mut particles = HashMap::new();

        for y in -1..=1 {
            for x in -1..=1 {
                let value = x as f64;
                let particle = Particle::new(x, y, params);
                particles.insert(particle, value);
            }
        }

        let map = ParticleMap::new(params, particles);
        let method = InterpolationMethod::IDW(IDWStrategy::default_from_params(&params));

        let test_cases = [
            // 1. Forward + ToDown
            GradTestCase {
                strategy: GradStrategy {
                    delta: 0.5,
                    difference_type: GradDifferenceType::Forward,
                    direction_type: GradDirectionType::ToDown,
                    ..Default::default()
                },
                name: "Forward+ToDown",
                expected_angles: vec![std::f64::consts::PI],
                value_threshold: -0.9,
                greater_than_threshold: false,
                check_abs: false,
            },
            // 2. Forward + ToUp
            GradTestCase {
                strategy: GradStrategy {
                    delta: 0.5,
                    difference_type: GradDifferenceType::Forward,
                    direction_type: GradDirectionType::ToUp,
                    ..Default::default()
                },
                name: "Forward+ToUp",
                expected_angles: vec![0.0],
                value_threshold: 0.9,
                greater_than_threshold: true,
                check_abs: false,
            },
            // 3. Forward + Steepest
            GradTestCase {
                strategy: GradStrategy {
                    delta: 0.5,
                    difference_type: GradDifferenceType::Forward,
                    direction_type: GradDirectionType::Steepest,
                    ..Default::default()
                },
                name: "Forward+Steepest",
                expected_angles: vec![0.0, std::f64::consts::PI],
                value_threshold: 0.9,
                greater_than_threshold: true,
                check_abs: true,
            },
            // 4. Central + ToDown
            GradTestCase {
                strategy: GradStrategy {
                    delta: 1.0,
                    difference_type: GradDifferenceType::Central,
                    direction_type: GradDirectionType::ToDown,
                    ..Default::default()
                },
                name: "Central+ToDown",
                expected_angles: vec![std::f64::consts::PI],
                value_threshold: -0.9,
                greater_than_threshold: false,
                check_abs: false,
            },
            // 5. Central + ToUp
            GradTestCase {
                strategy: GradStrategy {
                    delta: 1.0,
                    difference_type: GradDifferenceType::Central,
                    direction_type: GradDirectionType::ToUp,
                    ..Default::default()
                },
                name: "Central+ToUp",
                expected_angles: vec![0.0],
                value_threshold: 0.9,
                greater_than_threshold: true,
                check_abs: false,
            },
            // 6. Central + Steepest
            GradTestCase {
                strategy: GradStrategy {
                    delta: 1.0,
                    difference_type: GradDifferenceType::Central,
                    direction_type: GradDirectionType::Steepest,
                    ..Default::default()
                },
                name: "Central+Steepest",
                expected_angles: vec![0.0, std::f64::consts::PI],
                value_threshold: 0.9,
                greater_than_threshold: true,
                check_abs: true,
            },
        ];

        // Run all test cases
        for test_case in &test_cases {
            let gradient = map
                .get_gradient(0.0, 0.0, &test_case.strategy, &method)
                .unwrap();
            let gradient_unit = map.as_unit_vector(gradient.angle);

            // Check if the angle matches any of the expected angles
            let mut angle_matched = false;
            for &expected_angle in &test_case.expected_angles {
                let expected_unit = map.as_unit_vector(expected_angle);
                let unit_diff =
                    (gradient_unit.0 - expected_unit.0).hypot(gradient_unit.1 - expected_unit.1);

                if unit_diff < 0.1 {
                    angle_matched = true;
                    break;
                }
            }

            assert!(
                angle_matched,
                "{}: Expected angle near one of {:?}, got {:?}",
                test_case.name,
                test_case
                    .expected_angles
                    .iter()
                    .map(|&a| map.as_unit_vector(a))
                    .collect::<Vec<_>>(),
                gradient_unit
            );

            // Check gradient value
            let value_to_check = if test_case.check_abs {
                gradient.value.abs()
            } else {
                gradient.value
            };

            if test_case.greater_than_threshold {
                assert!(
                    value_to_check > test_case.value_threshold,
                    "{}: Expected gradient value greater than {}, got {}",
                    test_case.name,
                    test_case.value_threshold,
                    value_to_check
                );
            } else {
                assert!(
                    value_to_check < test_case.value_threshold,
                    "{}: Expected gradient value less than {}, got {}",
                    test_case.name,
                    test_case.value_threshold,
                    value_to_check
                );
            }
        }
    }

    #[test]
    fn test_get_gradient_no_interpolation() {
        let params = ParticleParameters::new(0.0, 0.0, 1.0, 0).unwrap();
        let map: ParticleMap<f64> = ParticleMap::new(params, HashMap::new()); // Empty map

        let strategy = GradStrategy {
            delta: 0.1,
            difference_type: GradDifferenceType::Forward,
            direction_type: GradDirectionType::ToDown,
            iteration: 2,
            sample_num: 8,
        };

        let gradient = map.get_gradient(0.0, 0.0, &strategy, &InterpolationMethod::Nearest);
        assert!(gradient.is_none(), "Expected None for empty map");
    }
}
