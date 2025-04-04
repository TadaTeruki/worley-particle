use worley_particle::{
    map::{algorithm::interp::nni::NNIStrategy, lerp::InterpolationMethod, ParticleMap},
    Particle, ParticleParameters,
};

fn main() {
    let params = ParticleParameters::new(0.9, 0.9, 0.5, 0).unwrap();
    let cells = Particle::from_inside_radius(0.0, 0.0, params, 5.0);
    let values = cells
        .iter()
        .map(|cell| (cell.hash_u64() % 10) as f64 * 0.1)
        .collect::<Vec<_>>();

    let map = ParticleMap::new(params, cells.into_iter().zip(values).collect());

    let nni_strategy = NNIStrategy::new_prebuild(&map);

    save_rasterised_map(
        &map,
        "data/output/map-rasterise-nni.png",
        300,
        300,
        nni_strategy.clone(),
    );
}

fn save_rasterised_map(
    map: &ParticleMap<f64>,
    filename: &str,
    image_width: usize,
    image_height: usize,
    nni_strategy: NNIStrategy,
) {
    let raster = map.rasterise(
        image_width,
        image_height,
        map.corners(),
        InterpolationMethod::NNI(nni_strategy),
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
