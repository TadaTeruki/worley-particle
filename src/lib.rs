/// Linear Grid ID.
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
            Self(x - 1, y),
            Self(x + 1, y),
            Self(x - 1, y + 1),
            Self(x, y + 1),
            Self(x + 1, y + 1),
        ]
    }
}

/// Non-Linear Grid ID.
/// x, y, and the degree of deformation (in \[0.0,1.0\]);
pub struct NLGridId(i64, i64, f64);

impl NLGridId {
    pub fn from(x: f64, y: f64, deformation: f64) -> Self {
        let lid = GridId::from(x, y);
        let around_lids = lid.get_around();

        let mut best_nlid = Self(lid.0, lid.1, deformation);
        let anchor = best_nlid.anchor_point();
        let mut best_sqdist = (anchor.0 as f64 - x).powi(2) + (anchor.1 as f64 - y).powi(2);
        for around_lid in around_lids {
            let around_nlid = Self(around_lid.0, around_lid.1, deformation);
            let anchor = around_nlid.anchor_point();
            let sqdist = (anchor.0 as f64 - x).powi(2) + (anchor.1 as f64 - y).powi(2);
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

    pub fn anchor_point(&self) -> (f64, f64) {
        let rel_core = anchor_point_from_hash(self.hash(), self.2);
        (self.0 as f64 + rel_core.0, self.1 as f64 + rel_core.1)
    }
}

fn xorshift64(x: u64) -> u64 {
    let mut x = x;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    x
}

fn hash_2d(x: u64, y: u64) -> u64 {
    let mut hash = x;
    hash = xorshift64(hash);
    hash = hash.wrapping_add(y as u64);
    hash = xorshift64(hash);
    hash
}

fn mod_f(a: f64, b: f64) -> f64 {
    a - b * (a / b).floor()
}

fn anchor_point_from_hash(hash: u64, diameter: f64) -> (f64, f64) {
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
                let point = anchor_point_from_hash(hash_2d(ix, iy), 1.0);
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
