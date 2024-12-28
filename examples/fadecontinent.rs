use std::collections::HashMap;

use nonlinear_grid::NLGridId;
struct FadePolygon {
    points: Vec<(f64, f64)>,
    center: (f64, f64),
}

fn grad(t: f64) -> f64 {
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

fn triangle_area(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1)).abs() / 2.0
}

fn point_in_triangle(p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    let s = triangle_area(a, b, c);
    let s1 = triangle_area(p, b, c);
    let s2 = triangle_area(a, p, c);
    let s3 = triangle_area(a, b, p);
    (s - s1 - s2 - s3).abs() < 1e-6
}

impl FadePolygon {
    pub fn new(points: Vec<(f64, f64)>) -> Self {
        let center = (
            points.iter().map(|p| p.0).sum::<f64>() / points.len() as f64,
            points.iter().map(|p| p.1).sum::<f64>() / points.len() as f64,
        );
        Self { points, center }
    }

    pub fn fade(&self, point: (f64, f64)) -> f64 {
        let rel_point = (point.0 - self.center.0, point.1 - self.center.1);
        let rel_poly = self
            .points
            .iter()
            .map(|p| (p.0 - self.center.0, p.1 - self.center.1))
            .collect::<Vec<(f64, f64)>>();
        let rel_center = (0.0, 0.0);

        // search triangle
        let mut triangle = (0, 0);
        for i in 0..self.points.len() {
            let j = (i + 1) % self.points.len();
            let a = rel_poly[i];
            let b = rel_poly[j];
            let c = rel_center;
            if point_in_triangle(rel_point, a, b, c) {
                triangle = (i, j);
                break;
            }
        }
        if triangle == (0, 0) {
            return 0.5;
        }

        let eu = rel_poly[triangle.1];
        let ev = rel_poly[triangle.0];
        let p = rel_point;
        let d = eu.0 * ev.1 - ev.0 * eu.1;

        let u = (-ev.0 * p.1 + ev.1 * p.0) / d;
        let v = (eu.0 * p.1 - eu.1 * p.0) / d;

        grad(1.0 - u - v)
    }
}

fn main() {
    let (min_x, max_x) = (-1.0, 1.0);
    let (min_y, max_y) = (-1.0, 1.0);
    let image_width = 500;
    let image_height = 500;

    let mut image_buf = image::RgbImage::new(image_width, image_height);

    let mut polygon_cache: HashMap<NLGridId, FadePolygon> = HashMap::new();

    let scale = 3.0;

    for iy in 0..image_height {
        for ix in 0..image_width {
            let x = min_x + (max_x - min_x) * ix as f64 / image_width as f64 * scale;
            let y = min_y + (max_y - min_y) * iy as f64 / image_height as f64 * scale;

            let nlgrid_id = NLGridId::from(x, y, 0.8);

            let fade = if let Some(polygon) = polygon_cache.get(&nlgrid_id) {
                polygon.fade((x, y))
            } else {
                let polygon = FadePolygon::new(nlgrid_id.calculate_voronoi().corners);
                let fade = polygon.fade((x, y));
                polygon_cache.insert(nlgrid_id, polygon);
                fade
            };
            let f = (fade * 255.0) as u8;
            image_buf.put_pixel(ix, iy, image::Rgb([f, f, f]));
        }
    }

    image_buf.save("out/image.png").unwrap();
}
