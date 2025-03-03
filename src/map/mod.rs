use std::{collections::HashMap, error::Error, fmt::Debug, io::Write};

use contour::ContourBuilder;
use flatbuffers::FlatBufferBuilder;
use particlemap_generated::particlemap::items::{
    finish_fbs_particle_map_buffer, root_as_fbs_particle_map, FbsParticle, FbsParticleArgs,
    FbsParticleMap, FbsParticleMapArgs, FbsParticleParameters, FbsParticleParametersArgs,
};
//use particlemap_generated::particlemap::items::root_as_particle_map;

use crate::{Particle, ParticleParameters};

#[path = "particlemap_generated.rs"]
#[allow(dead_code, unused_imports)]
mod particlemap_generated;

pub mod network;

static RW_SPLITTER: &str = "|";

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

pub trait ParticleMapAttributeRW: ParticleMapAttribute {
    fn from_strs(s: &[&str]) -> Result<Self, Box<dyn Error>>;
    fn to_strings(&self) -> Vec<String>;
    fn len_strs() -> usize;
}

impl ParticleMapAttributeRW for f64 {
    fn from_strs(s: &[&str]) -> Result<Self, Box<dyn Error>> {
        if let &[value] = s {
            Ok(value.parse()?)
        } else {
            Err("Expected one value".into())
        }
    }

    fn to_strings(&self) -> Vec<String> {
        vec![self.to_string()]
    }

    fn len_strs() -> usize {
        1
    }
}

impl ParticleMapAttributeRW for String {
    fn from_strs(s: &[&str]) -> Result<Self, Box<dyn Error>> {
        Ok(s.join(","))
    }

    fn to_strings(&self) -> Vec<String> {
        vec![self.clone()]
    }

    fn len_strs() -> usize {
        1
    }
}

impl ParticleMapAttributeRW for ParticleParameters {
    fn from_strs(s: &[&str]) -> Result<Self, Box<dyn Error>> {
        if let &[seed, min_randomness, max_randomness, scale] = s {
            Ok(Self {
                seed: seed.parse()?,
                min_randomness: min_randomness.parse()?,
                max_randomness: max_randomness.parse()?,
                scale: scale.parse()?,
            })
        } else {
            Err("Expected four values".into())
        }
    }

    fn to_strings(&self) -> Vec<String> {
        vec![
            self.seed.to_string(),
            self.min_randomness.to_string(),
            self.max_randomness.to_string(),
            self.scale.to_string(),
        ]
    }

    fn len_strs() -> usize {
        4
    }
}

impl ParticleMapAttributeRW for Particle {
    fn from_strs(s: &[&str]) -> Result<Self, Box<dyn Error>> {
        let s_coord = &s[0..2];
        let s_params = &s[2..];

        Ok(Self::new(
            s_coord[0].parse()?,
            s_coord[1].parse()?,
            ParticleParameters::from_strs(s_params)?,
        ))
    }

    fn to_strings(&self) -> Vec<String> {
        vec![self.grid.0.to_string(), self.grid.1.to_string()]
            .into_iter()
            .chain(self.params.to_strings().into_iter())
            .collect()
    }

    fn len_strs() -> usize {
        2 + ParticleParameters::len_strs()
    }
}

impl ParticleMapAttributeRW for () {
    fn from_strs(_: &[&str]) -> Result<Self, Box<dyn Error>> {
        Ok(())
    }

    fn to_strings(&self) -> Vec<String> {
        vec![]
    }

    fn len_strs() -> usize {
        0
    }
}

impl<T: ParticleMapAttributeRW> ParticleMap<T> {
    pub fn read_from_file(file_path: &str) -> Result<Self, Box<dyn Error>> {
        let data = std::fs::read(file_path)?;
        Self::read_from_bytes(data)
    }

    pub fn read_from_bytes(bytes: Vec<u8>) -> Result<Self, Box<dyn Error>> {
        let map = root_as_fbs_particle_map(&bytes)?;
        let fbs_params = map.params().ok_or("Missing parameters")?;
        let params = ParticleParameters {
            seed: fbs_params.seed(),
            min_randomness: fbs_params.min_randomness(),
            max_randomness: fbs_params.max_randomness(),
            scale: fbs_params.scale(),
        };

        let mut particles = HashMap::new();

        if let Some(fbs_particles) = map.particles() {
            for fbs_particle in fbs_particles.iter() {
                let value = if let Some(value_str) = fbs_particle.value() {
                    T::from_strs(&value_str.split(RW_SPLITTER).collect::<Vec<_>>())?
                } else {
                    return Err("Missing particle value".into());
                };
                particles.insert(
                    Particle::new(fbs_particle.x(), fbs_particle.y(), params),
                    value,
                );
            }
        }

        Ok(Self::new(params, particles))
    }

    pub fn write_to_bytes(&self) -> Result<Vec<u8>, Box<dyn Error>> {
        let mut builder = FlatBufferBuilder::new();

        // Create parameters
        let fbs_params = FbsParticleParameters::create(
            &mut builder,
            &FbsParticleParametersArgs {
                seed: self.params.seed,
                min_randomness: self.params.min_randomness,
                max_randomness: self.params.max_randomness,
                scale: self.params.scale,
            },
        );

        // Create particles
        let mut fbs_particles = Vec::new();
        for (particle, value) in &self.particles {
            let value_str = builder.create_string(&value.to_strings().join(RW_SPLITTER));
            let fbs_particle = FbsParticle::create(
                &mut builder,
                &FbsParticleArgs {
                    x: particle.grid.0,
                    y: particle.grid.1,
                    value: Some(value_str),
                },
            );
            fbs_particles.push(fbs_particle);
        }

        let fbs_particles_vec = builder.create_vector(&fbs_particles);

        // Create particle map
        let fbs_map = FbsParticleMap::create(
            &mut builder,
            &FbsParticleMapArgs {
                params: Some(fbs_params),
                particles: Some(fbs_particles_vec),
            },
        );

        finish_fbs_particle_map_buffer(&mut builder, fbs_map);

        Ok(builder.finished_data().to_vec())
    }

    pub fn write_to_file(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let data = self.write_to_bytes()?;
        let mut file = std::fs::File::create(file_path)?;
        file.write_all(&data)?;

        Ok(())
    }
}

pub trait ParticleMapAttributeLerp: ParticleMapAttribute + Default {
    fn lerp(&self, other: &Self, t: f64) -> Self;
}

impl ParticleMapAttributeLerp for f64 {
    fn lerp(&self, other: &Self, t: f64) -> Self {
        self + (other - self) * t
    }
}

impl ParticleMapAttributeLerp for () {
    fn lerp(&self, _: &Self, _: f64) -> Self {
        ()
    }
}

#[derive(Debug, Clone, PartialEq, Copy)]
pub struct IDWStrategy {
    pub sample_min_distance: f64,
    pub sample_max_distance: f64,
    pub weight_power: f64,
    pub smooth_power: Option<f64>,
}

impl IDWStrategy {
    pub fn default_from_params(params: &ParticleParameters) -> Self {
        Self {
            sample_min_distance: params.scale * 1e-6,
            sample_max_distance: params.scale * 1.415,
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

pub enum InterpolationMethod {
    Nearest,
    IDW(IDWStrategy),
    IDWSeparated(IDWStrategy),
}

impl<T: ParticleMapAttributeLerp> ParticleMap<T> {
    fn calculate_idw_weight(
        x: f64,
        y: f64,
        site: (f64, f64),
        sample_max_distance: f64,
        weight_power: f64,
        smooth_power: Option<f64>,
    ) -> Option<f64> {
        let distance = ((x - site.0).powi(2) + (y - site.1).powi(2)).sqrt();
        if distance > sample_max_distance {
            return None;
        }
        let weight = if let Some(smooth_power) = smooth_power {
            (1.0 - (distance / sample_max_distance).powf(smooth_power))
                / distance.powf(weight_power)
        } else {
            distance.powf(-weight_power)
        };

        Some(weight)
    }

    pub fn get_interpolated(&self, x: f64, y: f64, method: &InterpolationMethod) -> Option<T> {
        match method {
            InterpolationMethod::Nearest => {
                let particle = Particle::from(x, y, self.params);
                self.particles.get(&particle).cloned()
            }
            InterpolationMethod::IDW(strategy) => {
                let mut total_value: Option<T> = None;
                let mut tmp_weight = 0.0;
                let inside: Vec<Particle> =
                    Particle::from_inside_square(x, y, self.params, strategy.sample_max_distance);
                let particles_iter = inside.iter().filter_map(|particle| {
                    self.particles.get(particle).map(|value| (particle, value))
                });
                for (particle, value) in particles_iter {
                    let weight = Self::calculate_idw_weight(
                        x,
                        y,
                        particle.site(),
                        strategy.sample_max_distance,
                        strategy.weight_power,
                        strategy.smooth_power,
                    );
                    if let Some(weight) = weight {
                        tmp_weight += weight;
                        if let Some(total_value) = total_value.as_mut() {
                            *total_value = total_value.lerp(value, weight / tmp_weight);
                        } else {
                            total_value = Some(value.clone());
                        }
                    }
                }

                total_value
            }
            InterpolationMethod::IDWSeparated(strategy) => {
                let inside =
                    Particle::from_inside_square(x, y, self.params, strategy.sample_max_distance);
                let particles_iter = inside.iter().filter_map(|particle| {
                    self.particles.get(particle).map(|value| (particle, value))
                });
                let weights_iter = particles_iter.filter_map(|(particle, _)| {
                    Self::calculate_idw_weight(
                        x,
                        y,
                        particle.site(),
                        strategy.sample_max_distance,
                        strategy.weight_power,
                        strategy.smooth_power,
                    )
                    .map(|weight| (particle, weight))
                });

                let self_particle = Particle::from(x, y, self.params);

                let weight_sum = weights_iter.clone().map(|(_, weight)| weight).sum::<f64>();
                let mut total_value: Option<T> = None;
                let mut tmp_weight = 0.0;
                for (particle, weight) in weights_iter {
                    let (value, weight) = if particle == &self_particle {
                        (
                            self.particles.get(particle).unwrap(),
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

    pub fn rasterise(
        &self,
        img_width: usize,
        img_height: usize,
        corners: ((f64, f64), (f64, f64)),
        interp_method: &InterpolationMethod,
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
                *item = self.get_interpolated(x, y, interp_method);
            }
        }

        raster
    }
}

pub struct Band {
    pub threshold: f64,
    pub polygons: Vec<Vec<(f64, f64)>>,
}

enum VectorizationType {
    Contour,
    Isoband,
}

impl<T: ParticleMapAttributeLerp + Into<f64>> ParticleMap<T> {
    fn vectorize(
        &self,
        vectorization_type: VectorizationType,
        corners: ((f64, f64), (f64, f64)),
        rasterise_scale: f64,
        thresholds: &[f64],
        interp_method: &InterpolationMethod,
        smooth: bool,
    ) -> Result<Vec<Band>, Box<dyn Error>> {
        let ((original_min_x, original_min_y), (original_max_x, original_max_y)) = corners;

        let expansion = rasterise_scale * self.params.scale;

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
            interp_method,
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

        let result_coord_to_domain = |x: f64, y: f64| {
            (
                (x / width as f64) * (domain_max_x - domain_min_x) + domain_min_x,
                (y / height as f64) * (domain_max_y - domain_min_y) + domain_min_y,
            )
        };

        match vectorization_type {
            VectorizationType::Contour => {
                let result =
                    ContourBuilder::new(width, height, smooth).contours(&raster_f64, thresholds);
                if let Err(e) = result {
                    return Err(e.to_string().into());
                }

                let bands = result
                    .unwrap()
                    .iter()
                    .map(|contour| {
                        let polygons = contour
                            .geometry()
                            .iter()
                            .map(|polygon| {
                                let exterior = polygon.exterior();
                                exterior
                                    .0
                                    .iter()
                                    .map(|coord| result_coord_to_domain(coord.x, coord.y))
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>();

                        Band {
                            threshold: contour.threshold(),
                            polygons,
                        }
                    })
                    .collect::<Vec<_>>();

                Ok(bands)
            }
            VectorizationType::Isoband => {
                let result =
                    ContourBuilder::new(width, height, smooth).isobands(&raster_f64, thresholds);
                if let Err(e) = result {
                    return Err(e.to_string().into());
                }
                let bands = result
                    .unwrap()
                    .iter()
                    .map(|band| {
                        let polygons = band
                            .geometry()
                            .iter()
                            .map(|polygon| {
                                let exterior = polygon.exterior();
                                exterior
                                    .0
                                    .iter()
                                    .map(|coord| result_coord_to_domain(coord.x, coord.y))
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>();

                        Band {
                            threshold: band.min_v(),
                            polygons,
                        }
                    })
                    .collect::<Vec<_>>();

                Ok(bands)
            }
        }
    }

    pub fn contours(
        &self,
        corners: ((f64, f64), (f64, f64)),
        rasterise_scale: f64,
        thresholds: &[f64],
        interp_method: &InterpolationMethod,
        smooth: bool,
    ) -> Result<Vec<Band>, Box<dyn Error>> {
        self.vectorize(
            VectorizationType::Contour,
            corners,
            rasterise_scale,
            thresholds,
            interp_method,
            smooth,
        )
    }

    pub fn isobands(
        &self,
        corners: ((f64, f64), (f64, f64)),
        rasterise_scale: f64,
        thresholds: &[f64],
        interp_method: &InterpolationMethod,
        smooth: bool,
    ) -> Result<Vec<Band>, Box<dyn Error>> {
        self.vectorize(
            VectorizationType::Isoband,
            corners,
            rasterise_scale,
            thresholds,
            interp_method,
            smooth,
        )
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

    impl ParticleMapAttributeRW for TestAttribute {
        fn from_strs(s: &[&str]) -> Result<Self, Box<dyn Error>> {
            if let &[value, name] = s {
                Ok(Self {
                    value: value.parse()?,
                    name: name.to_string(),
                })
            } else {
                Err("Expected two values".into())
            }
        }

        fn to_strings(&self) -> Vec<String> {
            vec![self.value.to_string(), self.name.clone()]
        }

        fn len_strs() -> usize {
            2
        }
    }

    #[test]
    fn test_validate_read_and_write() {
        let params = ParticleParameters {
            seed: 0,
            min_randomness: 0.2,
            max_randomness: 0.5,
            scale: 0.1,
        };
        let particles = Particle::from_inside_radius(0.0, 0.0, params, 10.0);
        let values = particles
            .iter()
            .map(|particle| (particle.hash_u64() % 10) as f64)
            .collect::<Vec<_>>();
        let map = ParticleMap::new(params, particles.clone().into_iter().zip(values).collect());

        map.write_to_file("data/output/out.particlemap").unwrap();

        let map2 = ParticleMap::<f64>::read_from_file("data/output/out.particlemap").unwrap();

        assert!(map.is_same(&map2));

        let map = ParticleMap::new(
            params,
            particles
                .into_iter()
                .map(|particle| (particle, particle))
                .collect(),
        );

        map.write_to_file("data/output/out.particlemap").unwrap();

        let map2 = ParticleMap::<Particle>::read_from_file("data/output/out.particlemap").unwrap();

        assert!(map.is_same(&map2));
    }
}
