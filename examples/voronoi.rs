use std::collections::HashSet;

use tiny_skia::{Paint, PathBuilder, Pixmap, Stroke, Transform};
use worley_particle::{GenerationRules, Particle};

fn main() {
    let image_width = 500;
    let image_height = 500;
    let mut pixmap = Pixmap::new(image_width, image_height).unwrap();

    let (min_x, max_x) = (-1.0, 1.0);
    let (min_y, max_y) = (-1.0, 1.0);

    let mut cells = HashSet::new();

    for iy in 0..image_height {
        for ix in 0..image_width {
            let x = min_x + (max_x - min_x) * (ix as f64 / image_width as f64);
            let y = min_y + (max_y - min_y) * (iy as f64 / image_height as f64);

            let rules = GenerationRules::new(0.8, 0.8, 0.1, 0).unwrap();
            let ptc = Particle::from(x, y, rules);
            if cells.contains(&ptc) {
                continue;
            }
            cells.insert(ptc);
            let voronoi = ptc.calculate_voronoi().polygon;

            if voronoi.len() < 3 {
                continue;
            }

            let mut paint = Paint::default();
            let hash = ptc.hash_u64();
            let r = (hash % 256) as u8;
            let g = (hash / 256 % 256) as u8;
            let b = (hash / 65536 % 256) as u8;
            paint.set_color_rgba8(r, g, b, 255);

            let path = {
                let mut pb = PathBuilder::new();

                for i in 0..voronoi.len() {
                    let px = (voronoi[i].0 - min_x) / (max_x - min_x) * image_width as f64;
                    let py = (voronoi[i].1 - min_y) / (max_y - min_y) * image_height as f64;
                    if i == 0 {
                        pb.move_to(px as f32, py as f32);
                    } else {
                        pb.line_to(px as f32, py as f32);
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
    }
    pixmap.save_png("out/image.png").unwrap();
}
