use std::{
    cmp::Ordering,
    hash::{DefaultHasher, Hash, Hasher},
};

pub fn get_grids_around(x: i64, y: i64) -> [(i64, i64); 8] {
    [
        (x - 1, y - 1),
        (x, y - 1),
        (x + 1, y - 1),
        (x + 1, y),
        (x + 1, y + 1),
        (x, y + 1),
        (x - 1, y + 1),
        (x - 1, y),
    ]
}

pub fn get_grids_wide_around(x: i64, y: i64) -> [(i64, i64); 20] {
    [
        (x - 1, y - 2),
        (x, y - 2),
        (x + 1, y - 2),
        (x - 2, y - 1),
        (x - 1, y - 1),
        (x, y - 1),
        (x + 1, y - 1),
        (x + 2, y - 1),
        (x - 2, y),
        (x - 1, y),
        (x + 1, y),
        (x + 2, y),
        (x - 2, y + 1),
        (x - 1, y + 1),
        (x, y + 1),
        (x + 1, y + 1),
        (x + 2, y + 1),
        (x - 1, y + 2),
        (x, y + 2),
        (x + 1, y + 2),
    ]
}

pub fn square_distance(p1: &(f64, f64), p2: &(f64, f64)) -> f64 {
    (p1.0 - p2.0).powi(2) + (p1.1 - p2.1).powi(2)
}

pub fn circumcenter(triangle: &[&(f64, f64); 3]) -> (f64, f64) {
    let p1 = triangle[0];
    let p2 = triangle[1];
    let p3 = triangle[2];

    let d = 2.0 * (p1.0 * (p2.1 - p3.1) + p2.0 * (p3.1 - p1.1) + p3.0 * (p1.1 - p2.1));

    let ux = ((p1.0 * p1.0 + p1.1 * p1.1) * (p2.1 - p3.1)
        + (p2.0 * p2.0 + p2.1 * p2.1) * (p3.1 - p1.1)
        + (p3.0 * p3.0 + p3.1 * p3.1) * (p1.1 - p2.1))
        / d;
    let uy = ((p1.0 * p1.0 + p1.1 * p1.1) * (p3.0 - p2.0)
        + (p2.0 * p2.0 + p2.1 * p2.1) * (p1.0 - p3.0)
        + (p3.0 * p3.0 + p3.1 * p3.1) * (p2.0 - p1.0))
        / d;

    (ux, uy)
}

pub fn hash_2d(x: i64, y: i64, seed: u64) -> u64 {
    let mut hasher = DefaultHasher::new();
    x.hash(&mut hasher);
    y.hash(&mut hasher);
    seed.hash(&mut hasher);
    hasher.finish()
}

pub fn site_point_from_hash(hash: u64, min_diameter: f64, max_diameter: f64) -> (f64, f64) {
    let hash01 = (hash as f64 / u16::MAX as f64).fract();
    let theta = hash01 * 2.0 * std::f64::consts::PI;
    let radius = if min_diameter == max_diameter {
        min_diameter
    } else {
        min_diameter + (max_diameter - min_diameter) * hash01
    } / 2.0;
    (radius * theta.cos(), radius * theta.sin())
}

pub fn get_grid(x: f64, y: f64) -> (i64, i64) {
    (x.round() as i64, y.round() as i64)
}

pub fn arg_cmp(a: (f64, f64), b: (f64, f64)) -> Ordering {
    ((b.0, b.1) < (0., 0.))
        .cmp(&((a.0, a.1) < (0., 0.)))
        .then((a.0 * b.1).total_cmp(&(b.0 * a.1)))
}
