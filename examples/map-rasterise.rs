use worley_particle::{
    map::{
        lerp::{IDWStrategy, InterpolationMethod},
        ParticleMap,
    },
    Particle, ParticleParameters,
};

fn main() {
    let params = ParticleParameters::new(0.8, 0.8, 0.5, 0).unwrap();
    let cells = Particle::from_inside_radius(0.0, 0.0, params, 5.0);
    let values = cells
        .iter()
        .map(|cell| (cell.hash_u64() % 10) as f64 * 0.1)
        .collect::<Vec<_>>();

    let map = ParticleMap::new(params, cells.into_iter().zip(values).collect());

    save_rasterised_map(&map, "data/output/map-rasterise.png", 500, 500, &params);

    let rough_params = ParticleParameters::new(0.8, 0.8, 0.7, 0).unwrap();
    let rough_map = map
        .map_with_another_params_iter(
            InterpolationMethod::IDW(IDWStrategy::default_from_params(&params)),
            rough_params,
        )
        .collect::<ParticleMap<f64>>();

    save_rasterised_map(
        &rough_map,
        "data/output/map-rasterise-rough.png",
        500,
        500,
        &rough_params,
    );

    let fine_params = ParticleParameters::new(0.8, 0.8, 0.1, 0).unwrap();
    let fine_map = map
        .map_with_another_params_iter(
            InterpolationMethod::IDW(IDWStrategy::default_from_params(&params)),
            fine_params,
        )
        .collect::<ParticleMap<f64>>();

    save_rasterised_map(
        &fine_map,
        "data/output/map-rasterise-fine.png",
        500,
        500,
        &fine_params,
    );
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
