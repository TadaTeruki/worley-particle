use nonlinear_chunk::NLGridId;

fn main() {
    let (min_x, max_x) = (-1.0, 1.0);
    let (min_y, max_y) = (-1.0, 1.0);
    let image_width = 500;
    let image_height = 500;

    let mut image_buf = image::RgbImage::new(image_width, image_height);

    let scale = 3.0;

    for iy in 0..image_height {
        for ix in 0..image_width {
            let x = min_x + (max_x - min_x) * ix as f64 / image_width as f64 * scale;
            let y = min_y + (max_y - min_y) * iy as f64 / image_height as f64 * scale;

            let nlgrid_id = NLGridId::from(x, y, 0.8);
            let hash = nlgrid_id.hash();
            let r = (hash % 256) as u8;
            let g = (hash / 256 % 256) as u8;
            let b = (hash / 65536 % 256) as u8;

            image_buf.put_pixel(ix, iy, image::Rgb([r, g, b]));
        }
    }

    image_buf.save("out/image.png").unwrap();
}
