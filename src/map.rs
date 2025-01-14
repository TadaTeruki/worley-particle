use std::{collections::HashMap, error::Error, fmt::Debug};

use crate::{WorleyCell, WorleyParameters};

#[derive(Debug, Clone)]
pub struct WorleyMap<T: WorleyMapAttribute> {
    parameters: WorleyParameters,
    particles: HashMap<WorleyCell, T>,
}

pub trait WorleyMapAttribute: Debug + Clone + PartialEq {
    fn from_str(s: &[&str]) -> Result<Self, Box<dyn Error>>;
    fn to_string(&self) -> String;
}

pub trait WorleyMapAttributeLerp: WorleyMapAttribute {
    fn lerp(&self, other: &Self, t: f64) -> Self;
}

#[derive(Debug, Clone, PartialEq, Copy)]
pub struct IDWStrategy {
    pub sample_min_distance: f64,
    pub sample_max_distance: f64,
    pub weight_power: f64,
    pub smooth_power: Option<f64>,
}

impl IDWStrategy {
    pub fn default_from_parameters(parameters: &WorleyParameters) -> Self {
        Self {
            sample_min_distance: parameters.scale * 1e-6,
            sample_max_distance: parameters.scale * 1.414,
            weight_power: 1.5,
            smooth_power: Some(1.0),
        }
    }
}

impl Default for IDWStrategy {
    fn default() -> Self {
        Self {
            sample_min_distance: 1e-6,
            sample_max_distance: f64::INFINITY,
            weight_power: 1.0,
            smooth_power: Some(1.0),
        }
    }
}

pub enum RasterizationMethod {
    Nearest,
    IDW(IDWStrategy),
}

impl WorleyMapAttribute for f64 {
    fn from_str(s: &[&str]) -> Result<Self, Box<dyn Error>> {
        if let &[value] = s {
            Ok(value.parse()?)
        } else {
            Err("Expected one value".into())
        }
    }

    fn to_string(&self) -> String {
        ToString::to_string(self)
    }
}

impl WorleyMapAttributeLerp for f64 {
    fn lerp(&self, other: &Self, t: f64) -> Self {
        self + (other - self) * t
    }
}

impl WorleyMapAttribute for String {
    fn from_str(s: &[&str]) -> Result<Self, Box<dyn Error>> {
        Ok(s.join(","))
    }

    fn to_string(&self) -> String {
        self.clone()
    }
}

impl<T: WorleyMapAttribute> WorleyMap<T> {
    pub fn new(parameters: WorleyParameters, particles: HashMap<WorleyCell, T>) -> Self {
        Self {
            parameters,
            particles,
        }
    }

    pub fn read_from_file(file_path: &str) -> Result<Self, Box<dyn Error>> {
        let data = std::fs::read_to_string(file_path)?
            .lines()
            .map(|line| line.to_string())
            .collect();

        Self::read_from_lines(data)
    }

    pub fn read_from_lines(data: Vec<String>) -> Result<Self, Box<dyn Error>> {
        let mut parameters = WorleyParameters::default();

        let mut raw_particles: Vec<(f64, f64)> = Vec::new();
        let mut other_data: Vec<T> = Vec::new();

        for (i, line) in data.iter().enumerate() {
            let line = line.trim();
            if i == 0 {
                let iter = line.split(',');
                for value in iter {
                    let mut iter = value.split(':');
                    let key = iter.next().unwrap();
                    let value = iter.next().unwrap();
                    match key {
                        "seed" => parameters.seed = value.parse().unwrap(),
                        "min_randomness" => parameters.min_randomness = value.parse().unwrap(),
                        "max_randomness" => parameters.max_randomness = value.parse().unwrap(),
                        "scale" => parameters.scale = value.parse().unwrap(),
                        _ => panic!("Unknown key: {}", key),
                    }
                }
            } else {
                let mut iter = line.split('[');
                let mut iter = iter.nth(1).unwrap().split(']');
                let coords_str = iter.next().unwrap();
                let other_data_str = iter.next().unwrap();
                raw_particles.push({
                    let mut iter = coords_str.split(',');
                    let x = iter.next().unwrap().parse().unwrap();
                    let y = iter.next().unwrap().parse().unwrap();
                    (x, y)
                });
                other_data.push(T::from_str(&other_data_str.split(',').collect::<Vec<_>>())?);
            }
        }

        let particles = raw_particles
            .iter()
            .zip(other_data.iter())
            .map(|((x, y), other_data)| {
                let x = *x;
                let y = *y;
                let particle = WorleyCell::from(x, y, parameters);
                (particle, other_data.clone())
            })
            .collect::<HashMap<_, _>>();

        Ok(Self {
            parameters,
            particles,
        })
    }

    pub fn write_to_file(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut data = Vec::new();

        data.push(format!(
            "seed:{},min_randomness:{},max_randomness:{},scale:{}",
            self.particles.keys().next().unwrap().parameters.seed,
            self.particles
                .keys()
                .next()
                .unwrap()
                .parameters
                .min_randomness,
            self.particles
                .keys()
                .next()
                .unwrap()
                .parameters
                .max_randomness,
            self.particles.keys().next().unwrap().parameters.scale
        ));

        for (particle, other_data) in self.particles.iter() {
            let site = particle.site();
            data.push(format!("[{},{}]{}", site.0, site.1, other_data.to_string()));
        }

        std::fs::write(file_path, data.join("\n"))?;

        Ok(())
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

impl<T: WorleyMapAttribute + WorleyMapAttributeLerp> WorleyMap<T> {
    pub fn rasterise(
        &self,
        img_width: usize,
        img_height: usize,
        corners: ((f64, f64), (f64, f64)),
        method: RasterizationMethod,
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

                let particle = WorleyCell::from(x, y, self.parameters);

                let value = match method {
                    RasterizationMethod::Nearest => self.particles.get(&particle).cloned(),
                    RasterizationMethod::IDW(strategy) => {
                        let mut total_value: Option<T> = None;
                        let mut tmp_weight = 0.0;
                        let inside = WorleyCell::inside_square(
                            x,
                            y,
                            self.parameters,
                            strategy.sample_max_distance,
                        );
                        let particles_iter = inside
                            .iter()
                            .filter_map(|cell| self.particles.get(cell).map(|value| (cell, value)));
                        for (other_particle, value) in particles_iter {
                            let other_site = other_particle.site();
                            let distance =
                                ((x - other_site.0).powi(2) + (y - other_site.1).powi(2)).sqrt();
                            if distance > strategy.sample_max_distance {
                                continue;
                            }
                            if distance < strategy.sample_min_distance {
                                total_value = Some(value.clone());
                                break;
                            }
                            let weight = if let Some(smooth_power) = strategy.smooth_power {
                                (1.0 - (distance / strategy.sample_max_distance).powf(smooth_power))
                                    / distance.powf(strategy.weight_power)
                            } else {
                                distance.powf(-strategy.weight_power)
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
                };
                *item = value;
            }
        }

        raster
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Clone, PartialEq)]
    struct TestAttribute {
        value: f64,
        name: String,
    }

    impl WorleyMapAttribute for TestAttribute {
        fn from_str(s: &[&str]) -> Result<Self, Box<dyn Error>> {
            if let &[value, name] = s {
                Ok(Self {
                    value: value.parse()?,
                    name: name.to_string(),
                })
            } else {
                Err("Expected two values".into())
            }
        }

        fn to_string(&self) -> String {
            format!("{},{}", self.value, self.name)
        }
    }

    #[test]
    fn test_read_from_file() {
        let map = WorleyMap::<TestAttribute>::read_from_file("data/sample.worleymap").unwrap();

        assert_eq!(map.particles.len(), 3);
    }

    #[test]
    fn test_validate_read_and_write() {
        let parameters = WorleyParameters {
            seed: 0,
            min_randomness: 0.2,
            max_randomness: 0.5,
            scale: 0.1,
        };
        let cells = WorleyCell::inside_radius(0.0, 0.0, parameters, 10.0);
        let values = cells
            .iter()
            .map(|cell| (cell.hash_u64() % 10) as f64)
            .collect::<Vec<_>>();
        let map = WorleyMap::new(parameters, cells.into_iter().zip(values).collect());

        map.write_to_file("data/output/out.worleymap").unwrap();

        let map2 = WorleyMap::<f64>::read_from_file("data/test_cache/out.worleymap").unwrap();

        assert!(map.is_same(&map2));
    }
}
