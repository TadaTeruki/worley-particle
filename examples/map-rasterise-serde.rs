use std::{fs::File, path::Path};

use serde::{de::DeserializeOwned, Serialize};
use worley_particle::{
    map::{algorithm::interp::idw::IDWStrategy, lerp::InterpolationMethod, ParticleMap},
    Particle, ParticleParameters,
};

pub fn load_or_create<P: AsRef<Path>, T: Serialize + DeserializeOwned>(
    path: P,
    create_fn: impl FnOnce() -> T,
) -> Result<T, Box<dyn std::error::Error>> {
    let path = path.as_ref();

    if path.exists() {
        let file = File::open(path)?;
        let data: T = serde_json::from_reader(file)?;
        Ok(data)
    } else {
        let data = create_fn();
        let file = File::create(path)?;
        serde_json::to_writer(file, &data)?;
        Ok(data)
    }
}

fn create_map(params: ParticleParameters) -> ParticleMap<f64> {
    let cells = Particle::from_inside_radius(0.0, 0.0, params, 5.0);
    let values = cells
        .iter()
        .map(|cell| (cell.hash_u64() % 10) as f64 * 0.1)
        .collect::<Vec<_>>();

    ParticleMap::new(params, cells.into_iter().zip(values).collect())
}

#[cfg(feature = "serde")]
fn main() {
    let params = ParticleParameters::new(0.8, 0.8, 0.5, 0).unwrap();

    let map = load_or_create("data/output/map.json", || create_map(params.clone())).unwrap();

    save_rasterised_map(&map, "data/output/map-rasterise.png", 500, 500, &params);
}

fn save_rasterised_map(
    map: &ParticleMap<f64>,
    filename: &str,
    image_width: usize,
    image_height: usize,
    params: &ParticleParameters,
) {
    let raster = map.rasterise(
        image_width,
        image_height,
        map.corners(),
        InterpolationMethod::IDW(IDWStrategy::default_from_params(params)),
    );

    let mut image_buf = image::RgbImage::new(image_width as u32, image_height as u32);

    for (iy, row) in raster.iter().enumerate() {
        for (ix, value) in row.iter().enumerate() {
            if let Some(value) = value {
                let v = (value * 255.0) as u8;
                image_buf.put_pixel(ix as u32, iy as u32, image::Rgb([v, v, v]));
            } else {
                image_buf.put_pixel(ix as u32, iy as u32, image::Rgb([0, 0, 0]));
            }
        }
    }

    image_buf.save(filename).unwrap();
}
