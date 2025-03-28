use std::{collections::HashMap, error::Error, io::Write};

use flatbuffers::FlatBufferBuilder;
use particlemap_generated::particlemap::items::{
    finish_fbs_particle_map_buffer, root_as_fbs_particle_map, FbsParticle, FbsParticleArgs,
    FbsParticleMap, FbsParticleMapArgs, FbsParticleParameters, FbsParticleParametersArgs,
};

use crate::{Particle, ParticleParameters};

use super::{ParticleMap, ParticleMapAttribute};

#[path = "particlemap_generated.rs"]
#[allow(dead_code, unused_imports)]
#[allow(clippy::all)]
mod particlemap_generated;

static RW_SPLITTER: &str = "|";

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
            .chain(self.params.to_strings())
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
