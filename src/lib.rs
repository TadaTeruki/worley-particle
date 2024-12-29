use std::{
    cmp::Ordering,
    hash::{DefaultHasher, Hash, Hasher},
};

fn get_grids_around(x: i64, y: i64) -> [(i64, i64); 8] {
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

fn get_grids_wide_around(x: i64, y: i64) -> [(i64, i64); 20] {
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

fn square_distance(p1: &(f64, f64), p2: &(f64, f64)) -> f64 {
    (p1.0 - p2.0).powi(2) + (p1.1 - p2.1).powi(2)
}

fn circumcenter(triangle: &[&(f64, f64); 3]) -> (f64, f64) {
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

fn mod_f(a: f64, b: f64) -> f64 {
    a - b * (a / b).floor()
}

fn hash_2d(x: i64, y: i64) -> u64 {
    let mut hasher = DefaultHasher::new();
    x.hash(&mut hasher);
    y.hash(&mut hasher);
    hasher.finish()
}

fn site_point_from_hash(hash: u64, diameter: f64) -> (f64, f64) {
    let hash01 = mod_f(hash as f64 / u16::MAX as f64, 1.0);
    let theta = hash01 * 2.0 * std::f64::consts::PI;
    let radius = diameter / 2.0;
    (radius * theta.cos(), radius * theta.sin())
}

fn get_grid(x: f64, y: f64) -> (i64, i64) {
    (x.round() as i64, y.round() as i64)
}

fn arg_cmp(a: (f64, f64), b: (f64, f64)) -> Ordering {
    ((b.0, b.1) < (0., 0.))
        .cmp(&((a.0, a.1) < (0., 0.)))
        .then((a.0 * b.1).total_cmp(&(b.0 * a.1)))
}

/// A cell of worley noise.
/// The feature point (x, y) and the degree of deformation (in \[0.0,1.0\]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorleyCell(i64, i64, f64);

impl Hash for WorleyCell {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
        self.1.hash(state);
        self.2.to_bits().hash(state);
    }
}

impl Eq for WorleyCell {}

pub struct WorleyPolygon {
    pub polygon: Vec<(f64, f64)>,
    pub neighbors: Vec<WorleyCell>,
}

impl WorleyCell {
    pub fn new(x: i64, y: i64, deformation: f64) -> Self {
        Self(x, y, deformation)
    }

    pub fn from(x: f64, y: f64, deformation: f64) -> Self {
        let (ix, iy) = get_grid(x, y);
        let wc = Self(ix, iy, deformation);
        let mut best_wc = wc;
        let mut best_sqdist = square_distance(&(x, y), &wc.site());

        for (nx, ny) in get_grids_around(ix, iy).iter() {
            let wc = Self(*nx, *ny, deformation);
            let sqdist = square_distance(&(x, y), &wc.site());
            if sqdist < best_sqdist {
                best_wc = wc;
                best_sqdist = sqdist;
            }
        }

        best_wc
    }

    pub fn hash_u64(&self) -> u64 {
        hash_2d(self.0, self.1)
    }

    pub fn site(&self) -> (f64, f64) {
        let rel_core = site_point_from_hash(hash_2d(self.0, self.1), self.2);
        (self.0 as f64 + rel_core.0, self.1 as f64 + rel_core.1)
    }

    pub fn inside_radius(x: f64, y: f64, deformation: f64, radius: f64) -> Vec<Self> {
        let rad_ceil = radius.ceil() as i64;
        let mut surroundings = Vec::new();
        let (ix, iy) = get_grid(x, y);

        for dy in -rad_ceil..=rad_ceil {
            for dx in -rad_ceil..=rad_ceil {
                let wc = Self(ix + dx, iy + dy, deformation);
                if square_distance(&(x, y), &wc.site()) < radius.powi(2) {
                    surroundings.push(wc);
                }
            }
        }
        surroundings
    }

    /// Get the voronoi cell.
    pub fn calculate_voronoi(&self) -> WorleyPolygon {
        let site = self.site();
        let (ix, iy) = (self.0, self.1);

        // little bit complicated logic for large deformation
        let wide = self.2 > 0.5;

        let (mut around_sites, mut around_grids) = if wide {
            let around_grids = get_grids_wide_around(ix, iy);
            let around_sites = around_grids
                .iter()
                .map(|&x| Self(x.0, x.1, self.2).site())
                .collect::<Vec<(f64, f64)>>();

            let mut idx = (0..around_sites.len())
                .zip(around_sites.iter().map(|&x| (x.0 - site.0, x.1 - site.1)))
                .collect::<Vec<(usize, (f64, f64))>>();
            idx.sort_by(|a, b| arg_cmp(a.1, b.1));

            let around_sites = idx
                .iter()
                .map(|&(i, _)| around_sites[i])
                .collect::<Vec<(f64, f64)>>();
            let around_grids = idx
                .iter()
                .map(|&(i, _)| around_grids[i])
                .collect::<Vec<(i64, i64)>>();

            (around_sites, around_grids)
        } else {
            let around_grids = get_grids_around(ix, iy);

            let mut around_sites = vec![(0.0, 0.0); 8];
            for i in 0..around_grids.len() {
                around_sites[i] = Self(around_grids[i].0, around_grids[i].1, self.2).site();
            }
            (around_sites, around_grids.to_vec())
        };

        for _ in 0..65536 {
            let mut centers = vec![(0.0, 0.0); around_sites.len()];
            for i in 0..around_sites.len() {
                let ia = (i + 1) % around_sites.len();
                centers[i] = circumcenter(&[&around_sites[i], &around_sites[ia], &site]);
            }

            let mut applied = vec![false; around_sites.len()];
            for i in 0..around_sites.len() {
                let ib = (i + around_sites.len() - 1) % around_sites.len();
                let ia = (i + 1) % around_sites.len();

                let center_a = centers[i];
                let radius_a = square_distance(&center_a, &around_sites[i]);
                if square_distance(&center_a, &around_sites[ib]) < radius_a {
                    continue;
                }

                let center_b = centers[ib];
                let radius_b = square_distance(&center_b, &around_sites[i]);
                if square_distance(&center_b, &around_sites[ia]) < radius_b {
                    continue;
                }

                applied[i] = true;
            }

            (around_sites, around_grids) = around_sites
                .iter()
                .enumerate()
                .filter(|(i, _)| applied[*i])
                .map(|(i, &x)| (x, around_grids[i]))
                .unzip();

            if applied.iter().all(|&x| x) {
                break;
            }

            if !wide {
                break;
            }
        }

        let mut centers = Vec::new();

        for i in 0..around_sites.len() {
            let center = circumcenter(&[
                &around_sites[i],
                &around_sites[(i + 1) % around_sites.len()],
                &site,
            ]);
            centers.push(center);
        }

        WorleyPolygon {
            polygon: centers,
            neighbors: around_grids
                .iter()
                .map(|&x| Self(x.0, x.1, self.2))
                .collect(),
        }
    }
}

#[cfg(test)]
mod test {
    use rand::{rngs::StdRng, Rng, SeedableRng};

    use super::*;

    #[test]
    fn test_arg_cmp() {
        let points = {
            let mut rng = StdRng::from_seed([0; 32]);
            (0..100)
                .map(|_| (rng.gen_range(-100.0..100.0), rng.gen_range(-100.0..100.0)))
                .collect::<Vec<(f64, f64)>>()
        };
        let mut sorted_det = points.clone();
        let mut sorted_atan = points.clone();
        sorted_det.sort_by(|a, b| arg_cmp(*a, *b));
        sorted_atan.sort_by(|a, b| {
            let a = a.0.atan2(a.1);
            let b = b.0.atan2(b.1);
            a.total_cmp(&b)
        });

        // search same points
        let idx_atan = {
            let mut idx_atan = points.len();
            for i in 0..points.len() {
                if sorted_atan[i] == sorted_det[0] {
                    idx_atan = i;
                    break;
                }
            }
            idx_atan
        };

        assert!(idx_atan < points.len());

        for i in 0..points.len() {
            assert_eq!(sorted_det[i], sorted_atan[(i + idx_atan) % points.len()]);
        }
    }

    #[test]
    fn test_core_points() {
        let n = 1000;
        let mut xsum = 0.0;
        let mut ysum = 0.0;
        for iy in 0..n {
            for ix in 0..n {
                let point = site_point_from_hash(hash_2d(ix, iy), 1.0);
                xsum += point.0;
                ysum += point.1;
            }
        }

        let (xave, yave) = (xsum / (n * n) as f64, ysum / (n * n) as f64);
        let r = (xave.powi(2) + yave.powi(2)).sqrt();

        assert!(xave < 1e-3);
        assert!(yave < 1e-3);
        assert!(r < 1e-3);
    }

    #[test]
    fn test_calculate_voronoi() {
        let deform = 1.0 - 1e-9;
        let wc = WorleyCell::from(0.0, 0.0, deform);
        let voronoi = wc.calculate_voronoi();
        assert_eq!(voronoi.polygon.len(), voronoi.neighbors.len());

        for i in 0..voronoi.neighbors.len() {
            let a = voronoi.neighbors[i].site();
            let b = voronoi.neighbors[(i + 1) % voronoi.neighbors.len()].site();
            let c = wc.site();
            let center = circumcenter(&[&a, &b, &c]);
            assert!(square_distance(&center, &voronoi.polygon[i]) < 1e-9);
        }
    }

    #[test]
    fn test_inside_radius() {
        (0..100).for_each(|seed| {
            let mut rng = StdRng::from_seed([seed; 32]);
            let point: (f64, f64) = (rng.gen_range(-100.0..100.0), rng.gen_range(-100.0..100.0));
            let deform: f64 = rng.gen_range(0.0..1.0);
            let radius: f64 = rng.gen_range(0.0..100.0);
            let surroundings = WorleyCell::inside_radius(point.0, point.1, deform, radius);
            let mut count = 0;
            for iy in -(radius.ceil() + 1.0) as i64 * 5..=(radius.ceil() + 1.0) as i64 * 5 {
                for ix in -(radius.ceil() + 1.0) as i64 * 5..=(radius.ceil() + 1.0) as i64 * 5 {
                    let wc = WorleyCell::new(point.0 as i64 + ix, point.1 as i64 + iy, deform);
                    if square_distance(&point, &wc.site()) < radius.powi(2) {
                        count += 1;
                    }
                }
            }
            assert_eq!(surroundings.len(), count);
        });
    }
}
