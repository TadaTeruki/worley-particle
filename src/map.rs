use std::{collections::HashMap, error::Error, fmt::Debug};

use crate::{WorleyCell, WorleyParameters};

pub struct WorleyMap<T: WorleyMapAttribute> {
    particles: HashMap<WorleyCell, T>,
}

pub trait WorleyMapAttribute: Debug + Clone {
    fn from_str(s: &[&str]) -> Result<Self, Box<dyn Error>>;
    fn to_string(&self) -> String;
    fn lerp(&self, other: &Self, t: f64) -> Self;
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

    fn lerp(&self, other: &Self, t: f64) -> Self {
        if t < 0.5 {
            self.clone()
        } else {
            other.clone()
        }
    }
}

impl<T: WorleyMapAttribute> WorleyMap<T> {
    /// Read data from file
    pub fn read_from_file(file_path: &str) -> Result<Self, Box<dyn Error>> {
        let data = std::fs::read_to_string(file_path)?
            .lines()
            .map(|line| line.to_string())
            .collect();

        Self::read_from_vec(data)
    }

    pub fn read_from_vec(data: Vec<String>) -> Result<Self, Box<dyn Error>> {
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

        Ok(Self { particles })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Clone)]
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

        fn lerp(&self, other: &Self, t: f64) -> Self {
            Self {
                value: self.value + (other.value - self.value) * t,
                name: if t < 0.5 {
                    self.name.clone()
                } else {
                    other.name.clone()
                },
            }
        }
    }

    #[test]
    fn test_read_from_file() {
        let map = WorleyMap::<TestAttribute>::read_from_file("data/in.worleymap").unwrap();

        assert_eq!(map.particles.len(), 3);
    }
}
