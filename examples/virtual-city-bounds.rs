use worley_particle::{
    map::{lerp::InterpolationMethod, ParticleMap},
    Particle, ParticleParameters,
};

fn main() {
    let mut params = ParticleParameters::new(0.5, 0.5, 1.0, 0).unwrap();
    let cells = Particle::from_inside_radius(0.0, 0.0, params, 10.0);
    let values = cells
        .iter()
        .map(|cell| (cell.hash_u64() % 10) as f64 * 0.1)
        .collect::<Vec<_>>();

    let mut map = ParticleMap::new(params, cells.into_iter().zip(values).collect());

    for _ in 0..3 {
        params.scale *= 0.5;
        map = map
            .map_with_another_params_iter(InterpolationMethod::Nearest, params)
            .collect::<ParticleMap<f64>>();
    }

    save_rasterised_map(&map, "data/output/virtual-city-bounds.png", 500, 500);
}

fn save_rasterised_map(
    map: &ParticleMap<f64>,
    filename: &str,
    image_width: usize,
    image_height: usize,
) {
    let raster = map.rasterise(
        image_width,
        image_height,
        map.corners(),
        InterpolationMethod::Nearest,
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
