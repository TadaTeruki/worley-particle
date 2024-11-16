fn main() {
    let (min_x, max_x) = (-1.0, 1.0);
    let (min_y, max_y) = (-1.0, 1.0);
    let image_width = 800;
    let image_height = 800;

    let mut image_buf = image::RgbImage::new(image_width, image_height);

    for iy in 0..image_height {
        for ix in 0..image_width {
            let x = min_x + (max_x - min_x) * ix as f64 / image_width as f64;
            let y = min_y + (max_y - min_y) * iy as f64 / image_height as f64;

            let color = if x * x + y * y < 1.0 {
                image::Rgb([0, 0, 0])
            } else {
                image::Rgb([255, 255, 255])
            };

            image_buf.put_pixel(ix, iy, color);
        }
    }

    image_buf.save("out/image.png").unwrap();
}
