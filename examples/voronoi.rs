use nonlinear_chunk::NLGridId;
use tiny_skia::{Paint, PathBuilder, Pixmap, Stroke, Transform};

fn main() {
    let image_width = 500;
    let image_height = 500;
    let mut pixmap = Pixmap::new(image_width, image_height).unwrap();

    let n = 20;

    for iy in 0..n {
        for ix in 0..n {
            let nlgrid_id = NLGridId::from(ix as f64, iy as f64, 0.8);
            let voronoi = nlgrid_id.get_cell();

            if voronoi.len() < 3 {
                continue;
            }

            let mut paint = Paint::default();
            let hash = nlgrid_id.hash();
            let r = (hash % 256) as u8;
            let g = (hash / 256 % 256) as u8;
            let b = (hash / 65536 % 256) as u8;
            paint.set_color_rgba8(r, g, b, 255);

            let path = {
                let mut pb = PathBuilder::new();

                for i in 0..voronoi.len() {
                    let x = (voronoi[i].0) * image_width as f64 / (n - 1) as f64;
                    let y = (voronoi[i].1) * image_height as f64 / (n - 1) as f64;
                    if i == 0 {
                        pb.move_to(x as f32, y as f32);
                    } else {
                        pb.line_to(x as f32, y as f32);
                    }
                }
                pb.close();
                pb.finish().unwrap()
            };

            let mut stroke = Stroke::default();
            stroke.width = 2.0;

            pixmap.stroke_path(&path, &paint, &stroke, Transform::identity(), None);
        }
    }
    pixmap.save_png("out/image.png").unwrap();
}
