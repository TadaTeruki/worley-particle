use std::hash::{DefaultHasher, Hash, Hasher};

/// Linear Grid ID.
#[derive(Debug, Clone, Copy, PartialEq)]
struct GridId(i64, i64);

impl Eq for GridId {}

impl GridId {
    fn from(x: f64, y: f64) -> Self {
        Self(x.round() as i64, y.round() as i64)
    }

    fn get_around(&self) -> [Self; 8] {
        let (x, y) = (self.0, self.1);
        [
            Self(x - 1, y - 1),
            Self(x, y - 1),
            Self(x + 1, y - 1),
            Self(x + 1, y),
            Self(x + 1, y + 1),
            Self(x, y + 1),
            Self(x - 1, y + 1),
            Self(x - 1, y),
        ]
    }
    fn get_around_wide(&self) -> [Self; 20] {
        let (x, y) = (self.0, self.1);
        [
            Self(x - 1, y - 2),
            Self(x, y - 2),
            Self(x + 1, y - 2),
            Self(x - 2, y - 1),
            Self(x - 1, y - 1),
            Self(x, y - 1),
            Self(x + 1, y - 1),
            Self(x + 2, y - 1),
            Self(x - 2, y),
            Self(x - 1, y),
            Self(x + 1, y),
            Self(x + 2, y),
            Self(x - 2, y + 1),
            Self(x - 1, y + 1),
            Self(x, y + 1),
            Self(x + 1, y + 1),
            Self(x + 2, y + 1),
            Self(x - 1, y + 2),
            Self(x, y + 2),
            Self(x + 1, y + 2),
        ]
    }
}

/// Non-Linear Grid ID.
/// x, y, and the degree of deformation (in \[0.0,1.0\]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NLGridId(i64, i64, f64);

impl Hash for NLGridId {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
        self.1.hash(state);
        self.2.to_bits().hash(state);
    }
}

impl Eq for NLGridId {}

pub struct NLGridVoronoiCell {
    pub corners: Vec<(f64, f64)>,
    pub neighbors: Vec<NLGridId>,
}

impl NLGridId {
    pub fn from(x: f64, y: f64, deformation: f64) -> Self {
        let lid = GridId::from(x, y);
        let around_lids = lid.get_around();

        let mut best_nlid = Self(lid.0, lid.1, deformation);
        let site = best_nlid.site();
        let mut best_sqdist = square_distance(&(x, y), &site);
        for around_lid in around_lids {
            let around_nlid = Self(around_lid.0, around_lid.1, deformation);
            let site = around_nlid.site();
            let sqdist = square_distance(&(x, y), &site);
            if sqdist < best_sqdist {
                best_nlid = around_nlid;
                best_sqdist = sqdist;
            }
        }
        best_nlid
    }

    pub fn hash(&self) -> u64 {
        hash_2d(self.0, self.1)
    }

    pub fn site(&self) -> (f64, f64) {
        let rel_core = site_point_from_hash(hash_2d(self.0, self.1), self.2);
        (self.0 as f64 + rel_core.0, self.1 as f64 + rel_core.1)
    }

    /// Get the voronoi cell.
    pub fn calculate_voronoi(&self) -> NLGridVoronoiCell {
        let site = self.site();
        let lid = GridId::from(self.0 as f64, self.1 as f64);

        // little bit complicated logic for large deformation
        let wide = self.2 > 0.5;

        let (mut around_sites, mut around_lids) = if wide {
            let around_lids = lid.get_around_wide();

            let mut around_sites: Vec<(f64, f64)> = vec![(0.0, 0.0); 20];
            for i in 0..around_lids.len() {
                around_sites[i] = Self(around_lids[i].0, around_lids[i].1, self.2).site();
            }

            // sort by its phase
            let mut idx = (0..around_sites.len()).collect::<Vec<usize>>();
            idx.sort_by(|a, b| {
                let a_phase = (around_sites[*a].0 - site.0).atan2(around_sites[*a].1 - site.1);
                let b_phase = (around_sites[*b].0 - site.0).atan2(around_sites[*b].1 - site.1);
                a_phase.total_cmp(&b_phase)
            });

            let around_sites = idx
                .iter()
                .map(|&i| around_sites[i])
                .collect::<Vec<(f64, f64)>>();
            let around_lids = idx.iter().map(|&i| around_lids[i]).collect::<Vec<GridId>>();

            (around_sites, around_lids)
        } else {
            let around_lids = lid.get_around();

            let mut around_sites = vec![(0.0, 0.0); 8];
            for i in 0..around_lids.len() {
                around_sites[i] = Self(around_lids[i].0, around_lids[i].1, self.2).site();
            }
            (around_sites, around_lids.to_vec())
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

            (around_sites, around_lids) = around_sites
                .iter()
                .enumerate()
                .filter(|(i, _)| applied[*i])
                .map(|(i, &x)| (x, around_lids[i]))
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

        NLGridVoronoiCell {
            corners: centers,
            neighbors: around_lids
                .iter()
                .map(|&x| Self(x.0, x.1, self.2))
                .collect(),
        }
    }
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

#[cfg(test)]
mod test {
    use super::*;

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
        let nlid = NLGridId::from(0.0, 0.0, deform);
        let voronoi = nlid.calculate_voronoi();
        assert_eq!(voronoi.corners.len(), voronoi.neighbors.len());

        for i in 0..voronoi.neighbors.len() {
            let a = voronoi.neighbors[i].site();
            let b = voronoi.neighbors[(i + 1) % voronoi.neighbors.len()].site();
            let c = nlid.site();
            let center = circumcenter(&[&a, &b, &c]);
            assert!(square_distance(&center, &voronoi.corners[i]) < 1e-9);
        }
    }
}
