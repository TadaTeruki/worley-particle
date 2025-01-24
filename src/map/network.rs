use std::collections::{BinaryHeap, HashMap};

use crate::{Particle, ParticleParameters};

use super::{ParticleMap, ParticleMapAttribute};

pub struct ParticleNetwork {
    params: ParticleParameters,
    network: HashMap<(i64, i64), Vec<(i64, i64)>>,
}

impl ParticleNetwork {
    pub fn new<T: ParticleMapAttribute>(map: &ParticleMap<T>) -> Self {
        let mut network = HashMap::new();
        map.iter().for_each(|(particle, _)| {
            let neighbors = particle.calculate_voronoi().neighbors;
            let key = (particle.grid.0, particle.grid.1);
            for neighbor in neighbors {
                let neighbor_key = (neighbor.grid.0, neighbor.grid.1);
                network
                    .entry(neighbor_key)
                    .or_insert_with(Vec::new)
                    .push(key);
            }
        });
        Self {
            params: *map.params(),
            network,
        }
    }

    pub fn network_into_hashmap(&self) -> HashMap<Particle, Vec<Particle>> {
        self.network
            .clone()
            .into_iter()
            .map(|(k, v)| {
                let k_particle = Particle::new(k.0, k.1, self.params);
                let mut v_particles = Vec::new();
                for neighbor in v {
                    let neighbor_particle = Particle::new(neighbor.0, neighbor.1, self.params);
                    v_particles.push(neighbor_particle);
                }
                (k_particle, v_particles)
            })
            .collect::<HashMap<Particle, Vec<Particle>>>()
    }

    pub fn is_connected(&self, a: Particle, b: Particle) -> bool {
        if let Some(neighbors) = self.network.get(&a.grid) {
            neighbors.contains(&b.grid)
        } else {
            false
        }
    }

    pub fn neighbors(&self, a: Particle) -> Vec<Particle> {
        if let Some(neighbors) = self.network.get(&a.grid) {
            neighbors
                .iter()
                .map(|n| Particle::new(n.0, n.1, self.params))
                .collect()
        } else {
            Vec::new()
        }
    }

    pub fn create_path(
        &self,
        start: Particle,
        mut finish_condition: impl FnMut(Particle) -> bool,
        mut evaluation: impl FnMut(Particle, Particle) -> f64,
    ) -> Vec<Particle> {
        let mut heap = BinaryHeap::new();
        heap.push(StateForPathfinding {
            particle: start,
            cost: 0.0,
        });
        let mut parent = HashMap::new();
        parent.insert(start, start);

        loop {
            let current_state = heap.pop().unwrap();

            let neighbors = self.neighbors(current_state.particle);
            let mut next = None;
            let mut next_evaluation = std::f64::MAX;
            for neighbor in neighbors {
                if parent.contains_key(&neighbor) {
                    continue;
                }
                let evaluation = evaluation(current_state.particle, neighbor);
                if evaluation > next_evaluation {
                    next = Some(neighbor);
                    next_evaluation = evaluation;
                }
            }
            if let Some(next) = next {
                heap.push(StateForPathfinding {
                    particle: next,
                    cost: next_evaluation,
                });
                parent.insert(next, current_state.particle);
            } else {
                break;
            }
            if finish_condition(current_state.particle) {
                break;
            }
        }

        let mut path = Vec::new();
        let mut current = start;
        while current != parent[&current] {
            path.push(current);
            current = parent[&current];
        }
        path.push(start);
        path.reverse();
        path
    }
}

#[derive(Debug, Clone, PartialEq)]
struct StateForPathfinding {
    particle: Particle,
    cost: f64,
}

impl Eq for StateForPathfinding {}

impl PartialOrd for StateForPathfinding {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.cost.partial_cmp(&other.cost)
    }
}

impl Ord for StateForPathfinding {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cost.partial_cmp(&other.cost).unwrap()
    }
}
