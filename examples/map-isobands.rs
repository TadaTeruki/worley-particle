use tiny_skia::{Paint, PathBuilder, Pixmap, Stroke, Transform};
use worley_particle::{
    map::{IDWStrategy, IsobandResult, ParticleMap, RasteriseMethod},
    GenerationRules, Particle,
};

fn main() {
    let rules = GenerationRules::new(0.8, 0.8, 0.5, 0).unwrap();
    let radius = 5.0;
    let cells = Particle::from_inside_radius(2.0, 1.0, rules, radius);
    let values = cells
        .iter()
        .map(|cell| (cell.hash_u64() % 10) as f64 * 0.1)
        .collect::<Vec<_>>();

    let map = ParticleMap::new(rules, cells.into_iter().zip(values).collect());

    let thresholds = (0..=10).map(|i| i as f64 / 10.).collect::<Vec<_>>();

    let ((min_x, min_y), (max_x, max_y)) = map.corners();

    let isobands = map
        .isobands(
            ((min_x, min_y), (max_x, max_y)),
            30.0,
            &thresholds,
            RasteriseMethod::IDW(IDWStrategy::default_from_rules(&rules)),
            true,
        )
        .unwrap();

    let mut paint = Paint::default();

    let image_width = 500;
    let image_height = 500;

    let mut pixmap = Pixmap::new(image_width, image_height).unwrap();

    for IsobandResult {
        threshold,
        polygons,
    } in isobands
    {
        let grayscale = (threshold * 255.0) as u8;
        paint.set_color_rgba8(grayscale, grayscale, grayscale, 255);
        let path = {
            let mut pb = PathBuilder::new();
            for polygon in polygons {
                for i in 0..polygon.len() {
                    let px = (polygon[i].0 - min_x) / (max_x - min_x) * image_width as f64;
                    let py = (polygon[i].1 - min_y) / (max_y - min_y) * image_height as f64;
                    if i == 0 {
                        pb.move_to(px as f32, py as f32);
                    } else {
                        pb.line_to(px as f32, py as f32);
                    }
                }
            }
            pb.close();
            pb.finish()
        };

        if path.is_none() {
            continue;
        }

        let mut stroke = Stroke::default();
        stroke.width = 2.0;

        pixmap.stroke_path(&path.unwrap(), &paint, &stroke, Transform::identity(), None);
    }

    pixmap.save_png("data/output/map-isobands.png").unwrap();

    // --- Overlay the isobands map on top of the rasterised map ---

    let isoband_image = image::open("data/output/map-isobands.png")
        .unwrap()
        .to_rgb8();

    let raster = map.rasterise(
        image_width as usize,
        image_height as usize,
        map.corners(),
        RasteriseMethod::IDW(IDWStrategy::default_from_rules(&rules)),
    );

    let mut image_buf = image::RgbImage::new(image_width as u32, image_height as u32);

    for y in 0..image_height {
        for x in 0..image_width {
            let img_pixel = isoband_image.get_pixel(x as u32, y as u32).0;
            let img_pixel = image::Rgb([img_pixel[0], img_pixel[1], img_pixel[2]]);

            let raster_pixel = raster[y as usize][x as usize];
            if raster_pixel.is_none() {
                continue;
            }
            let raster_pixel = raster_pixel.unwrap();
            let raster_pixel = (raster_pixel * 255.0) as u8;
            let img_pixel = image::Rgb([
                (img_pixel[0] as f64 * 0.2 + raster_pixel as f64) as u8,
                (img_pixel[1] as f64 * 0.2 + raster_pixel as f64) as u8,
                (img_pixel[2] as f64 * 0.2 + raster_pixel as f64) as u8,
            ]);

            image_buf.put_pixel(x as u32, y as u32, img_pixel);
        }
    }

    image_buf
        .save("data/output/map-isobands-overlay.png")
        .unwrap();
}
