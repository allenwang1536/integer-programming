use std::{
    cmp::{Ordering, Reverse},
    collections::BinaryHeap,
    f64, usize,
};

use ordered_float::NotNan;

use crate::{
    ffi::{
        status::{INFEASIBLE, OPTIMAL},
        OrConstraintHandle, OrSolverHandle, OrVarHandle,
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
}

struct Fixing {
    var: OrVarHandle,
    val: f64,
}

struct Node {
    parent: NodeId,
    fixing: Option<Fixing>,
    lb: f64,
}

#[derive(Clone, Copy, Eq, PartialEq)]
struct SearchHead {
    node: NodeId,
    lb: NotNan<f64>,
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
        }
    }

    fn reset_vars(&mut self) {
        for v in &self.test_var {
            v.set_bounds(0.0, 1.0);
        }
    }

    fn insert_node(&mut self, node: Node) -> NodeId {
        let id = self.nodes.len() as NodeId;
        self.nodes.push(node);
        id
    }

    fn apply_fixings(&mut self, node: NodeId) {
        let mut cur = node;
        while cur != NO_PARENT {
            if let Some(Fixing { var, val }) = &self.nodes[cur].fixing {
                var.set_bounds(*val, *val);
            }
            cur = self.nodes[cur].parent;
        }
    }

    pub fn solve(&mut self) -> usize {
        let root_id = self.insert_node(Node {
            parent: NO_PARENT,
            fixing: None,
            lb: 0.0,
        });
        self.heads.push(SearchHead {
            node: root_id,
            lb: NotNan::new(0.0).unwrap(),
        });

        let mut ip_opt = f64::INFINITY;
        let mut soln = vec![0u8; self.instance.num_tests];

        while let Some(cur) = self.heads.pop() {
            self.reset_vars();
            self.apply_fixings(cur.node);

            // do solve
            match self.lp_solver.solve() {
                OPTIMAL => {
                    // if cur LB is worse than incumb, drop
                    let cur_lb = self.lp_solver.objective_value();
                    if cur_lb >= ip_opt {
                        continue;
                    }
                    // find non-integral var
                    if let Some(var) = self.find_non_integral() {
                        // branch on this var
                        let l = self.insert_node(Node {
                            parent: cur.node,
                            lb: cur_lb,
                            fixing: Some(Fixing {
                                var: var.clone(),
                                val: 0.0,
                            }),
                        });
                        let r = self.insert_node(Node {
                            parent: cur.node,
                            lb: cur_lb,
                            fixing: Some(Fixing { var, val: 1.0 }),
                        });
                        self.heads.push(SearchHead {
                            node: r,
                            lb: NotNan::new(cur_lb).unwrap(),
                        });
                        self.heads.push(SearchHead {
                            node: l,
                            lb: NotNan::new(cur_lb).unwrap(),
                        });
                    } else {
                        // this is ip feasible
                        println!("found some ip feasible with obj: {}", cur_lb);
                        if cur_lb < ip_opt {
                            ip_opt = cur_lb;
                            for (t, v) in self.test_var.iter().enumerate() {
                                soln[t] = v.solution_value() as u8;
                            }
                        }
                    }
                }
                INFEASIBLE => {
                    // nooo
                }
                _ => unreachable!(),
            }
        }

        ip_opt as usize
    }

    pub fn find_non_integral(&self) -> Option<OrVarHandle> {
        for t in &self.test_var {
            match t.solution_value() {
                0.0 | 1.0 => continue,
                _ => return Some(t.clone()),
            };
        }
        None
    }

    pub fn print_soln(&self) {
        for (t, v) in self.test_var.iter().enumerate() {
            println!("use_test[{}] = {}", t, v.solution_value())
        }
    }
}
