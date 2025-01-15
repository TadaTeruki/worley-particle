use std::{collections::HashMap, error::Error, fmt::Debug};

use contour_isobands::isobands;

use crate::{GenerationRules, Particle};

#[derive(Debug, Clone)]
pub struct ParticleMap<T: Debug + Clone + PartialEq> {
    rules: GenerationRules,
    particles: HashMap<Particle, T>,
}

pub trait ParticleMapRWAttribute: Debug + Clone + PartialEq {
    fn from_str(s: &[&str]) -> Result<Self, Box<dyn Error>>;
    fn to_string(&self) -> String;
}

pub trait ParticleMapAttributeLerp: ParticleMapRWAttribute {
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
    pub fn default_from_rules(rules: &GenerationRules) -> Self {
        Self {
            sample_min_distance: rules.scale * 1e-6,
            sample_max_distance: rules.scale * 1.415,
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
            weight_power: 1.5,
            smooth_power: Some(1.0),
        }
    }
}

pub enum RasteriseMethod {
    Nearest,
    IDW(IDWStrategy),
}

impl ParticleMapRWAttribute for f64 {
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

impl ParticleMapAttributeLerp for f64 {
    fn lerp(&self, other: &Self, t: f64) -> Self {
        self + (other - self) * t
    }
}

impl ParticleMapRWAttribute for String {
    fn from_str(s: &[&str]) -> Result<Self, Box<dyn Error>> {
        Ok(s.join(","))
    }

    fn to_string(&self) -> String {
        self.clone()
    }
}

impl<T: ParticleMapRWAttribute> ParticleMap<T> {
    pub fn read_from_file(file_path: &str) -> Result<Self, Box<dyn Error>> {
        let data = std::fs::read_to_string(file_path)?
            .lines()
            .map(|line| line.to_string())
            .collect();

        Self::read_from_lines(data)
    }

    pub fn read_from_lines(data: Vec<String>) -> Result<Self, Box<dyn Error>> {
        let mut rules = GenerationRules::default();

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
                        "seed" => rules.seed = value.parse().unwrap(),
                        "min_randomness" => rules.min_randomness = value.parse().unwrap(),
                        "max_randomness" => rules.max_randomness = value.parse().unwrap(),
                        "scale" => rules.scale = value.parse().unwrap(),
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
                let particle = Particle::from(x, y, rules);
                (particle, other_data.clone())
            })
            .collect::<HashMap<_, _>>();

        Ok(Self { rules, particles })
    }

    pub fn write_to_file(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let mut data = Vec::new();

        data.push(format!(
            "seed:{},min_randomness:{},max_randomness:{},scale:{}",
            self.particles.keys().next().unwrap().rules.seed,
            self.particles.keys().next().unwrap().rules.min_randomness,
            self.particles.keys().next().unwrap().rules.max_randomness,
            self.particles.keys().next().unwrap().rules.scale
        ));

        for (particle, other_data) in self.particles.iter() {
            let site = particle.site();
            data.push(format!("[{},{}]{}", site.0, site.1, other_data.to_string()));
        }

        std::fs::write(file_path, data.join("\n"))?;

        Ok(())
    }
}

impl<T: Debug + Clone + PartialEq> ParticleMap<T> {
    pub fn new(rules: GenerationRules, particles: HashMap<Particle, T>) -> Self {
        Self { rules, particles }
    }

    pub fn corners(&self) -> ((f64, f64), (f64, f64)) {
        let particle_sites = self
            .particles
            .keys()
            .map(|cell| cell.site())
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

impl<T: Debug + Clone + PartialEq + ParticleMapAttributeLerp> ParticleMap<T> {
    pub fn rasterise(
        &self,
        img_width: usize,
        img_height: usize,
        corners: ((f64, f64), (f64, f64)),
        method: RasteriseMethod,
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

                let particle = Particle::from(x, y, self.rules);

                let value = match method {
                    RasteriseMethod::Nearest => self.particles.get(&particle).cloned(),
                    RasteriseMethod::IDW(strategy) => {
                        let mut total_value: Option<T> = None;
                        let mut tmp_weight = 0.0;
                        let inside = Particle::from_inside_square(
                            x,
                            y,
                            self.rules,
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

pub struct IsobandResult {
    pub threshold: f64,
    pub polygons: Vec<Vec<(f64, f64)>>,
}

impl<T: Debug + Clone + PartialEq + ParticleMapAttributeLerp + Into<f64>> ParticleMap<T> {
    pub fn isobands(
        &self,
        corners: ((f64, f64), (f64, f64)),
        rasterise_scale: f64,
        thresholds: &[f64],
        rasterise_method: RasteriseMethod,
        multi_thread: bool,
    ) -> Result<Vec<IsobandResult>, Box<dyn Error>> {
        let ((original_min_x, original_min_y), (original_max_x, original_max_y)) = corners;

        let expansion = rasterise_scale * self.rules.scale;

        let raster_min_x = (original_min_x * expansion).floor() as i64;
        let raster_max_x = (original_max_x * expansion).ceil() as i64;
        let raster_min_y = (original_min_y * expansion).floor() as i64;
        let raster_max_y = (original_max_y * expansion).ceil() as i64;

        let domain_min_x = raster_min_x as f64 / expansion;
        let domain_max_x = raster_max_x as f64 / expansion;
        let domain_min_y = raster_min_y as f64 / expansion;
        let domain_max_y = raster_max_y as f64 / expansion;

        let width = (raster_max_x - raster_min_x) as usize;
        let height = (raster_max_y - raster_min_y) as usize;
        let raster = self.rasterise(
            width,
            height,
            ((domain_min_x, domain_min_y), (domain_max_x, domain_max_y)),
            rasterise_method,
        );

        let raster_f64 = raster
            .iter()
            .flat_map(|row| {
                row.iter()
                    .map(|value: &Option<T>| {
                        if let Some(value) = value {
                            value.clone().into()
                        } else {
                            f64::MIN
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        let res_raster = isobands(
            &raster_f64,
            thresholds,
            multi_thread,
            width,
            height,
            multi_thread,
        );
        if let Err(e) = res_raster {
            return Err(e.to_string().into());
        }

        let res_domain = res_raster
            .unwrap()
            .iter()
            .map(|band| {
                let threshold = band.1;
                let polygons = band
                    .0
                    .iter()
                    .map(|polygon| {
                        polygon
                            .iter()
                            .map(|point| {
                                let x = (point.x() / width as f64) * (domain_max_x - domain_min_x)
                                    + domain_min_x;
                                let y = (point.y() / height as f64) * (domain_max_y - domain_min_y)
                                    + domain_min_y;
                                (x, y)
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                IsobandResult {
                    threshold,
                    polygons,
                }
            })
            .collect::<Vec<_>>();

        Ok(res_domain)
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

    impl ParticleMapRWAttribute for TestAttribute {
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
        let map = ParticleMap::<TestAttribute>::read_from_file("data/sample.particlemap").unwrap();

        assert_eq!(map.particles.len(), 3);
    }

    #[test]
    fn test_validate_read_and_write() {
        let rules = GenerationRules {
            seed: 0,
            min_randomness: 0.2,
            max_randomness: 0.5,
            scale: 0.1,
        };
        let cells = Particle::from_inside_radius(0.0, 0.0, rules, 10.0);
        let values = cells
            .iter()
            .map(|cell| (cell.hash_u64() % 10) as f64)
            .collect::<Vec<_>>();
        let map = ParticleMap::new(rules, cells.into_iter().zip(values).collect());

        map.write_to_file("data/output/out.particlemap").unwrap();

        let map2 = ParticleMap::<f64>::read_from_file("data/output/out.particlemap").unwrap();

        assert!(map.is_same(&map2));
    }
}
