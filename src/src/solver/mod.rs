use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    f64,
    time::{Duration, Instant},
    usize,
};

use ordered_float::NotNan;

use crate::{
    ffi::{
        status::{INFEASIBLE, OPTIMAL},
        OrBasisHandle, OrConstraintHandle, OrSolverHandle, OrVarHandle,
    },
    instance::IPInstance,
};

type NodeId = usize;
const NO_PARENT: NodeId = usize::MAX;

pub struct BNCSolver {
    lp_solver: OrSolverHandle,
    test_var: Vec<OrVarHandle>,
    constraints: Vec<OrConstraintHandle>,
    nodes: Vec<Node>,
    heads: BinaryHeap<SearchHead>,
    pub instance: IPInstance,
    pub time: Duration,
}

struct Fixing {
    var: u32,
    val: u8,
}

struct Node {
    parent: NodeId,
    fixing: Option<Fixing>,
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum NodeStatus {
    Integral,
    Fractional { branch_var: u32 },
}

#[derive(Clone, Eq, PartialEq)]
struct SearchHead {
    node: NodeId,
    lb: NotNan<f64>,
    status: NodeStatus,
    basis: OrBasisHandle,
}

impl Ord for SearchHead {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .lb
            .cmp(&self.lb)
            .then_with(|| self.node.cmp(&other.node))
    }
}

impl PartialOrd for SearchHead {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl BNCSolver {
    pub fn new(instance: IPInstance) -> Self {
        let lp_solver = OrSolverHandle::new_scip();
        let mut test_var = Vec::with_capacity(instance.num_tests);
        let num_constaints = instance.num_diseases * (instance.num_diseases - 1) / 2;
        let mut constraints = Vec::with_capacity(num_constaints);

        for t in 0..instance.num_tests {
            let name = format!("t{}", t);
            test_var.push(lp_solver.make_num_var(0.0, 1.0, &name))
        }

        for d_l in 0..instance.num_diseases {
            for d_r in d_l + 1..instance.num_diseases {
                let c = lp_solver.make_constraint(1.0, f64::INFINITY);
                for (t, v) in test_var.iter().enumerate() {
                    let diff = instance.at(t, d_l) ^ instance.at(t, d_r);
                    lp_solver.constraint_set_coefficient(&c, v, diff as f64);
                }
                constraints.push(c)
            }
        }

        for (t, c) in test_var.iter().zip(instance.test_costs.iter()) {
            lp_solver.objective_set_coefficient(t, *c as f64);
        }
        lp_solver.objective_set_minimize();

        BNCSolver {
            lp_solver,
            test_var,
            constraints,
            nodes: Vec::new(),
            heads: BinaryHeap::new(),
            instance,
            time: Duration::new(0, 0),
        }
    }

    fn reset_vars(&mut self) {
        let t = Instant::now();
        for v in &self.test_var {
            v.set_bounds(0.0, 1.0);
        }
        self.time += t.elapsed();
    }

    fn insert_node(&mut self, node: Node) -> NodeId {
        let id = self.nodes.len();
        self.nodes.push(node);
        id
    }

    fn insert_child(&mut self, parent: NodeId, var: u32, val: u8) -> NodeId {
        let id = self.nodes.len();

        self.nodes.push(Node {
            parent,
            fixing: Some(Fixing { var, val }),
        });

        id
    }

    fn apply_fixings(&mut self, node: NodeId) {
        let mut cur = node;
        while cur != NO_PARENT {
            if let Some(Fixing { var, val }) = self.nodes[cur].fixing {
                self.test_var[var as usize].set_bounds(val as f64, val as f64);
            }
            cur = self.nodes[cur].parent;
        }
    }

    fn solve_node(&mut self, node: NodeId, basis: Option<&OrBasisHandle>) -> Option<SearchHead> {
        self.reset_vars();
        self.apply_fixings(node);

        if let Some(basis) = basis {
            self.lp_solver.restore_basis(basis);
        }

        match self.lp_solver.solve() {
            INFEASIBLE => None,
            OPTIMAL => {
                let lb = NotNan::new(self.lp_solver.objective_value()).unwrap();

                // find non-integral var
                let status = match self.find_non_integral() {
                    Some(var) => NodeStatus::Fractional { branch_var: var },
                    None => NodeStatus::Integral,
                };

                Some(SearchHead {
                    node,
                    lb,
                    status,
                    basis: self.lp_solver.save_basis(),
                })
            }
            _ => unreachable!(),
        }
    }

    pub fn solve(&mut self) -> usize {
        let root_id = self.insert_node(Node {
            parent: NO_PARENT,
            fixing: None,
        });

        self.lp_solver.solve(); // assuming there exists a solution here, wtv
        let Some(branch_var) = self.find_non_integral() else {
            return self.lp_solver.objective_value() as usize;
        };

        self.heads.push(SearchHead {
            node: root_id,
            lb: NotNan::new(self.lp_solver.objective_value()).unwrap(),
            status: NodeStatus::Fractional { branch_var },
            basis: self.lp_solver.save_basis(),
        });

        let mut ip_opt = NotNan::new(f64::INFINITY).unwrap();
        let mut soln = vec![0u8; self.instance.num_tests];

        // let mut i = 0;

        while let Some(cur) = self.heads.pop() {
            // if i % 10 == 0 {
            //     println!(
            //         "num nodes: {}, num_heads: {}",
            //         self.nodes.len(),
            //         self.heads.len()
            //     )
            // }
            // i += 1;

            // if cur LB is worse than incumb, drop
            if cur.lb >= ip_opt {
                break;
            }

            match cur.status {
                NodeStatus::Integral => {
                    // this is ip feasible
                    if cur.lb < ip_opt {
                        ip_opt = cur.lb;
                        for (t, v) in self.test_var.iter().enumerate() {
                            soln[t] = v.solution_value() as u8;
                        }
                    }
                }
                NodeStatus::Fractional { branch_var } => {
                    let l = self.insert_child(cur.node, branch_var, 0);
                    let r = self.insert_child(cur.node, branch_var, 1);
                    if let Some(l_head) = self.solve_node(l, Some(&cur.basis)) {
                        self.heads.push(l_head);
                    }
                    if let Some(r_head) = self.solve_node(r, Some(&cur.basis)) {
                        self.heads.push(r_head);
                    }
                }
            }
        }

        *ip_opt as usize
    }

    pub fn find_non_integral(&self) -> Option<u32> {
        let mut most_frac = None;
        for t in 0..self.instance.num_tests {
            match self.test_var[t].solution_value() {
                0.0 | 1.0 => continue,
                n => {
                    let cur_dist = (n - 0.5).abs();
                    match most_frac {
                        Some((_, dist)) => {
                            if cur_dist < dist {
                                most_frac = Some((t as u32, cur_dist))
                            }
                        }

                        None => most_frac = Some((t as u32, (n - 0.5).abs())),
                    }
                }
            };
        }

        most_frac.map(|v| v.0)
    }

    pub fn print_soln(&self) {
        for (t, v) in self.test_var.iter().enumerate() {
            println!("use_test[{}] = {}", t, v.solution_value())
        }
    }
}
