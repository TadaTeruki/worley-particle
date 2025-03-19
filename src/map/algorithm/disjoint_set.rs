use std::collections::HashMap;

use crate::Particle;

pub struct DisjointSet {
    parent: HashMap<Particle, Particle>,
    rank: HashMap<Particle, usize>,
}

impl DisjointSet {
    pub fn new(xs: &[Particle]) -> Self {
        let parent = xs.iter().map(|&x| (x, x)).collect::<HashMap<_, _>>();
        let rank = xs.iter().map(|&x| (x, 0)).collect::<HashMap<_, _>>();
        Self { parent, rank }
    }

    pub fn find(&mut self, x: Particle) -> Option<Particle> {
        if let Some(&parent) = self.parent.get(&x) {
            if parent == x {
                Some(x)
            } else {
                let root = self.find(parent)?;
                self.parent.insert(x, root);
                Some(root)
            }
        } else {
            None
        }
    }

    pub fn union(&mut self, x: Particle, y: Particle) -> bool {
        let x_root = match self.find(x) {
            Some(root) => root,
            None => return false,
        };

        let y_root = match self.find(y) {
            Some(root) => root,
            None => return false,
        };

        if x_root == y_root {
            return true;
        }

        match (self.rank.get(&x_root), self.rank.get(&y_root)) {
            (Some(&x_rank), Some(&y_rank)) => {
                if x_rank < y_rank {
                    self.parent.insert(x_root, y_root);
                } else if x_rank > y_rank {
                    self.parent.insert(y_root, x_root);
                } else {
                    self.parent.insert(y_root, x_root);
                    self.rank.insert(x_root, x_rank + 1);
                }
                true
            }
            _ => false,
        }
    }

    pub fn is_same_set(&mut self, x: Particle, y: Particle) -> bool {
        match (self.find(x), self.find(y)) {
            (Some(x_root), Some(y_root)) => x_root == y_root,
            _ => false,
        }
    }

    fn get_set_elements(&mut self, x: Particle) -> Vec<Particle> {
        let x_root = match self.find(x) {
            Some(root) => root,
            None => return Vec::new(),
        };

        let keys: Vec<Particle> = self.parent.keys().cloned().collect();
        keys.into_iter()
            .filter(|&k| self.find(k) == Some(x_root))
            .collect()
    }

    fn get_all_roots(&mut self) -> Vec<Particle> {
        let mut roots = Vec::new();
        let keys: Vec<Particle> = self.parent.keys().cloned().collect();
        for x in keys {
            if let Some(root) = self.find(x) {
                if !roots.contains(&root) {
                    roots.push(root);
                }
            }
        }
        roots
    }

    pub fn get_all_sets(&mut self) -> Vec<Vec<Particle>> {
        let roots = self.get_all_roots();

        roots
            .into_iter()
            .map(|root| self.get_set_elements(root))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_disjoint_set() {
        let params = crate::ParticleParameters::default();

        let p1 = Particle::new(0, 0, params);
        let p2 = Particle::new(1, 0, params);
        let p3 = Particle::new(0, 1, params);
        let p4 = Particle::new(1, 1, params);

        let mut ds = DisjointSet::new(&[p1, p2, p3, p4]);

        assert!(!ds.is_same_set(p1, p2));
        assert!(!ds.is_same_set(p1, p3));
        assert!(!ds.is_same_set(p1, p4));

        ds.union(p1, p2);
        assert!(ds.is_same_set(p1, p2));
        assert!(!ds.is_same_set(p1, p3));

        ds.union(p3, p4);
        assert!(ds.is_same_set(p3, p4));
        assert!(!ds.is_same_set(p1, p3));

        ds.union(p1, p3);
        assert!(ds.is_same_set(p1, p3));
        assert!(ds.is_same_set(p1, p4));
        assert!(ds.is_same_set(p2, p3));

        let sets = ds.get_all_sets();
        assert_eq!(sets.len(), 1);
        assert_eq!(sets[0].len(), 4);
    }
}
