use std::f64;

use crate::{
    ffi::{
        status::{INFEASIBLE, OPTIMAL},
        OrConstraintHandle, OrSolverHandle, OrVarHandle,
    },
    instance::IPInstance,
};

pub struct BNCSolver {
    lp_solver: OrSolverHandle,
    test_var: Vec<OrVarHandle>,
    constraints: Vec<OrConstraintHandle>,
    pub instance: IPInstance,
}

struct Fixing {
    var: OrVarHandle,
    val: f64,
}

enum Node {
    Solve(Option<Fixing>),
    Undo(OrVarHandle),
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
            instance,
        }
    }

    pub fn solve(&mut self) -> usize {
        let mut stack = vec![Node::Solve(None)];

        let mut ip_opt = f64::INFINITY;
        let mut soln = vec![0u8; self.instance.num_tests];

        while let Some(cur) = stack.pop() {
            match cur {
                Node::Solve(fixing) => {
                    // fix
                    if let Some(Fixing { var, val }) = fixing {
                        var.set_bounds(val, val);
                    }

                    // do solve
                    match self.lp_solver.solve() {
                        OPTIMAL => {
                            // if cur LB is worse than incumb, drop
                            let cur_obj = self.lp_solver.objective_value();
                            if cur_obj >= ip_opt {
                                continue;
                            }
                            // find non-integral var
                            if let Some(var) = self.find_non_integral() {
                                // branch on this var
                                stack.push(Node::Undo(var.clone()));
                                stack.push(Node::Solve(Some(Fixing {
                                    var: var.clone(),
                                    val: 1.0,
                                })));
                                stack.push(Node::Solve(Some(Fixing { var, val: 0.0 })));
                            } else {
                                // this is ip feasible
                                println!("found some ip feasible with obj: {}", cur_obj);
                                if cur_obj < ip_opt {
                                    ip_opt = cur_obj;
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

                Node::Undo(var) => {
                    var.set_bounds(0.0, 1.0);
                }
            }
        }

        // for (t, v) in soln.iter().enumerate() {
        //     println!("use_test[{}] = {}", t, v);
        // }

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
