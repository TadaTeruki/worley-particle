/// Linear Grid ID.
#[derive(Debug)]
struct GridId(i64, i64);

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
#[derive(Debug)]
pub struct NLGridId(i64, i64, f64);

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
        hash_2d(self.0 as u64, self.1 as u64)
    }

    pub fn site(&self) -> (f64, f64) {
        let rel_core = site_point_from_hash(self.hash(), self.2);
        (self.0 as f64 + rel_core.0, self.1 as f64 + rel_core.1)
    }

    /// Get the voronoi cell.
    pub fn get_cell(&self) -> Vec<(f64, f64)> {
        let site = self.site();
        let lid = GridId::from(self.0 as f64, self.1 as f64);

        // little bit complicated logic for large deformation
        let wide = self.2 > 0.5;

        let mut around_sites = if wide {
            let around_lids = lid.get_around_wide();
            
            let mut around_sites = vec![(0.0, 0.0); 20];
            for i in 0..around_lids.len() {
                around_sites[i] = Self(around_lids[i].0, around_lids[i].1, self.2).site();
            }
            // sort by its phase
            around_sites.sort_by(|a, b| {
                let a_phase = (a.0 - site.0).atan2(a.1 - site.1);
                let b_phase = (b.0 - site.0).atan2(b.1 - site.1);
                a_phase.total_cmp(&b_phase)
            });

            around_sites
        } else {
            let around_lids = lid.get_around();
            
            let mut around_sites = vec![(0.0, 0.0); 8];
            for i in 0..around_lids.len() {
                around_sites[i] = Self(around_lids[i].0, around_lids[i].1, self.2).site();
            }
            around_sites
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

            around_sites = around_sites
                .iter()
                .enumerate()
                .filter(|(i, _)| applied[*i])
                .map(|(_, &x)| x)
                .collect::<Vec<(f64, f64)>>();

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
        centers
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

fn xorshift64(x: u64) -> u64 {
    let mut x = x;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    x
}

fn hash_2d(x: u64, y: u64) -> u64 {
    let mut hash = x.wrapping_mul(y+17320508);
    hash = xorshift64(hash);
    hash = hash.wrapping_add(y);
    hash = xorshift64(hash);
    hash
}

fn mod_f(a: f64, b: f64) -> f64 {
    a - b * (a / b).floor()
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
}
