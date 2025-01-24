use std::{
    cmp::Ordering,
    hash::{DefaultHasher, Hash, Hasher},
};

use thiserror::Error;

#[cfg(feature = "particlemap")]
pub mod map;

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

fn hash_2d(x: i64, y: i64, seed: u64) -> u64 {
    let mut hasher = DefaultHasher::new();
    x.hash(&mut hasher);
    y.hash(&mut hasher);
    seed.hash(&mut hasher);
    hasher.finish()
}

fn site_point_from_hash(hash: u64, min_diameter: f64, max_diameter: f64) -> (f64, f64) {
    let hash01 = (hash as f64 / u16::MAX as f64).fract();
    let theta = hash01 * 2.0 * std::f64::consts::PI;
    let radius = if min_diameter == max_diameter {
        min_diameter
    } else {
        min_diameter + (max_diameter - min_diameter) * hash01
    } / 2.0;
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

#[derive(Debug, Error)]
pub enum GenerationRuleError {
    #[error("Invalid randomness (randomnesses are not in [0.0, 1.0], or min_randomness > max_randomness")]
    InvalidRandomness,

    #[error("Invalid scale (scale is not positive)")]
    InvalidScale,
}

/// Parameters for the generation of the Worley noise.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParticleParameters {
    /// Minimum randomness of the feature point (in [0.0, 1.0]).
    /// randomness <= 0.5 is recommended for better performance and stability.
    min_randomness: f64,
    /// Maximum randomness of the feature point (in [0.0, 1.0]).
    /// randomness <= 0.5 is recommended for better performance and stability.
    max_randomness: f64,
    /// Scale of the grid.
    scale: f64,
    /// The seed of the random number generator.
    seed: u64,
}

impl Hash for ParticleParameters {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.min_randomness.to_bits().hash(state);
        self.max_randomness.to_bits().hash(state);
        self.scale.to_bits().hash(state);
        self.seed.hash(state);
    }
}

impl Eq for ParticleParameters {}

impl Default for ParticleParameters {
    fn default() -> Self {
        Self {
            min_randomness: 0.5,
            max_randomness: 0.5,
            scale: 1.0,
            seed: 0,
        }
    }
}

impl ParticleParameters {
    pub fn new(
        min_randomness: f64,
        max_randomness: f64,
        scale: f64,
        seed: u64,
    ) -> Result<Self, GenerationRuleError> {
        Self::default()
            .set_randomness(min_randomness, max_randomness)?
            .set_scale(scale)?
            .set_seed(seed)
    }

    pub fn set_randomness(
        mut self,
        min_randomness: f64,
        max_randomness: f64,
    ) -> Result<Self, GenerationRuleError> {
        if min_randomness > max_randomness || min_randomness < 0.0 || max_randomness > 1.0 {
            return Err(GenerationRuleError::InvalidRandomness);
        }
        self.min_randomness = min_randomness;
        self.max_randomness = max_randomness;
        Ok(self)
    }

    pub fn set_scale(mut self, scale: f64) -> Result<Self, GenerationRuleError> {
        if scale <= 0.0 {
            return Err(GenerationRuleError::InvalidScale);
        }
        self.scale = scale;
        Ok(self)
    }

    pub fn set_seed(mut self, seed: u64) -> Result<Self, GenerationRuleError> {
        self.seed = seed;
        Ok(self)
    }

    pub fn get_randomness(&self) -> (f64, f64) {
        (self.min_randomness, self.max_randomness)
    }

    pub fn get_scale(&self) -> f64 {
        self.scale
    }
}

pub struct ParticlePolygon {
    pub polygon: Vec<(f64, f64)>,
    pub neighbors: Vec<Particle>,
}

impl ParticlePolygon {
    pub fn area(&self) -> f64 {
        let mut area = 0.0;
        for i in 0..self.polygon.len() {
            let a = self.polygon[i];
            let b = self.polygon[(i + 1) % self.polygon.len()];
            area += a.0 * b.1 - a.1 * b.0;
        }
        area.abs() / 2.0
    }
}

/// A particle in the Worley noise.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Particle {
    /// (x, y) coordinate of the grid.
    grid: (i64, i64),
    /// rule of the generation.
    params: ParticleParameters,
}

impl Hash for Particle {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.grid.hash(state);
        self.params.hash(state);
    }
}

impl Eq for Particle {}

impl Particle {
    pub fn new(grid_x: i64, grid_y: i64, params: ParticleParameters) -> Self {
        Self {
            grid: (grid_x, grid_y),
            params,
        }
    }

    pub fn from(x: f64, y: f64, params: ParticleParameters) -> Self {
        let (ix, iy) = get_grid(x / params.scale, y / params.scale);
        let ptc = Self::new(ix, iy, params);
        let mut best_ptc = ptc;
        let mut best_sqdist = square_distance(&(x, y), &ptc.site());

        for (nx, ny) in get_grids_around(ix, iy).iter() {
            let ptc = Self::new(*nx, *ny, params);
            let sqdist = square_distance(&(x, y), &ptc.site());
            if sqdist < best_sqdist {
                best_ptc = ptc;
                best_sqdist = sqdist;
            }
        }

        best_ptc
    }

    pub fn hash_u64(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }

    pub fn site(&self) -> (f64, f64) {
        let rel_core = site_point_from_hash(
            hash_2d(self.grid.0, self.grid.1, self.params.seed),
            self.params.min_randomness,
            self.params.max_randomness,
        );
        (
            (self.grid.0 as f64 + rel_core.0) * self.params.scale,
            (self.grid.1 as f64 + rel_core.1) * self.params.scale,
        )
    }

    pub fn from_inside_radius(
        x: f64,
        y: f64,
        params: ParticleParameters,
        radius: f64,
    ) -> Vec<Self> {
        let (corner_min_x, corner_min_y) =
            get_grid((x - radius) / params.scale, (y - radius) / params.scale);
        let (corner_max_x, corner_max_y) =
            get_grid((x + radius) / params.scale, (y + radius) / params.scale);

        let mut surroundings = Vec::new();
        for iy in corner_min_y..=corner_max_y {
            for ix in corner_min_x..=corner_max_x {
                let ptc = Self::new(ix, iy, params);
                if square_distance(&(x, y), &ptc.site()) < radius.powi(2) {
                    surroundings.push(ptc);
                }
            }
        }

        surroundings
    }

    pub fn from_inside_square(x: f64, y: f64, params: ParticleParameters, side: f64) -> Vec<Self> {
        let (corner_min_x, corner_min_y) =
            get_grid((x - side) / params.scale, (y - side) / params.scale);
        let (corner_max_x, corner_max_y) =
            get_grid((x + side) / params.scale, (y + side) / params.scale);

        let mut surroundings = Vec::new();
        for iy in corner_min_y..=corner_max_y {
            for ix in corner_min_x..=corner_max_x {
                surroundings.push(Self::new(ix, iy, params));
            }
        }

        surroundings
    }

    /// Get the voronoi cell.
    pub fn calculate_voronoi(&self) -> ParticlePolygon {
        let site = self.site();
        let (ix, iy) = (self.grid.0, self.grid.1);

        // little bit complicated logic for large max_randomness
        let wide = self.params.max_randomness > 0.5;

        let (mut around_sites, mut around_grids) = if wide {
            let around_grids = get_grids_wide_around(ix, iy);
            let around_sites = around_grids
                .iter()
                .map(|&x| Self::new(x.0, x.1, self.params).site())
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
                around_sites[i] =
                    Self::new(around_grids[i].0, around_grids[i].1, self.params).site();
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

        ParticlePolygon {
            polygon: centers,
            neighbors: around_grids
                .iter()
                .map(|&x| Self::new(x.0, x.1, self.params))
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
    // check if the core points are distributed uniformly
    fn test_core_points() {
        let n = 1000;
        let mut xsum = 0.0;
        let mut ysum = 0.0;
        for iy in 0..n {
            for ix in 0..n {
                let point = site_point_from_hash(hash_2d(ix, iy, 0), 1.0, 1.0);
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
        let randomness = 0.9;
        let params = ParticleParameters::new(randomness, randomness, 1.0, 0).unwrap();
        let ptc = Particle::from(0.0, 0.0, params);
        let voronoi = ptc.calculate_voronoi();
        assert_eq!(voronoi.polygon.len(), voronoi.neighbors.len());

        for i in 0..voronoi.neighbors.len() {
            let a = voronoi.neighbors[i].site();
            let b = voronoi.neighbors[(i + 1) % voronoi.neighbors.len()].site();
            let c = ptc.site();
            let center = circumcenter(&[&a, &b, &c]);
            assert!(square_distance(&center, &voronoi.polygon[i]) < 1e-9);
        }
    }

    #[test]
    fn test_inside_radius() {
        (0..50).for_each(|seed| {
            let mut rng = StdRng::from_seed([seed; 32]);
            let point: (f64, f64) = (rng.gen_range(-100.0..100.0), rng.gen_range(-100.0..100.0));
            let randomness: f64 = rng.gen_range(0.0..0.9);
            let params = ParticleParameters::new(randomness, randomness, 1.0, 0).unwrap();
            let radius: f64 = rng.gen_range(0.0..100.0);
            let mut count = 0;
            let surroundings = Particle::from_inside_radius(point.0, point.1, params, radius);
            for iy in -(radius.ceil() + 1.0) as i64 * 3..=(radius.ceil() + 1.0) as i64 * 3 {
                for ix in -(radius.ceil() + 1.0) as i64 * 3..=(radius.ceil() + 1.0) as i64 * 3 {
                    let ptc = Particle::new(point.0 as i64 + ix, point.1 as i64 + iy, params);
                    if square_distance(&point, &ptc.site()) < radius.powi(2) {
                        count += 1;
                    }
                }
            }
            assert_eq!(surroundings.len(), count);
        });
    }
}
