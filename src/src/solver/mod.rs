use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    f64,
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering as AtomicOrd},
        Condvar, Mutex, RwLock,
    },
    thread, usize,
};

use ordered_float::NotNan;

use crate::{
    ffi::{
        status::{INFEASIBLE, OPTIMAL},
        OrBasisHandle, OrSolverHandle, OrVarHandle,
    },
    instance::IPInstance,
};

unsafe impl Send for OrBasisHandle {}
unsafe impl Send for OrSolverHandle {}
unsafe impl Send for OrVarHandle {}

type NodeId = usize;
const NO_PARENT: NodeId = usize::MAX;

pub struct BNCSolver {
    pub instance: IPInstance,
}

#[derive(Clone, Copy)]
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

struct SharedState {
    nodes: RwLock<Vec<Node>>,
    heads: Mutex<BinaryHeap<SearchHead>>,
    heads_cv: Condvar,
    ip_opt: Mutex<NotNan<f64>>,
    active: AtomicUsize,
    done: AtomicBool,
}

struct WorkerSolver {
    lp_solver: OrSolverHandle,
    test_var: Vec<OrVarHandle>,
}

impl WorkerSolver {
    fn new(instance: &IPInstance) -> Self {
        let lp_solver = OrSolverHandle::new_scip();
        let mut test_var = Vec::with_capacity(instance.num_tests);

        for t in 0..instance.num_tests {
            let name = format!("t{}", t);
            test_var.push(lp_solver.make_num_var(0.0, 1.0, &name));
        }

        for d_l in 0..instance.num_diseases {
            for d_r in d_l + 1..instance.num_diseases {
                let c = lp_solver.make_constraint(1.0, f64::INFINITY);
                for (t, v) in test_var.iter().enumerate() {
                    let diff = instance.at(t, d_l) ^ instance.at(t, d_r);
                    lp_solver.constraint_set_coefficient(&c, v, diff as f64);
                }
            }
        }

        for (t, c) in test_var.iter().zip(instance.test_costs.iter()) {
            lp_solver.objective_set_coefficient(t, *c as f64);
        }
        lp_solver.objective_set_minimize();

        WorkerSolver {
            lp_solver,
            test_var,
        }
    }

    fn reset_and_apply(&self, fixings: &[(u32, u8)]) {
        for v in &self.test_var {
            v.set_bounds(0.0, 1.0);
        }
        for &(var, val) in fixings {
            self.test_var[var as usize].set_bounds(val as f64, val as f64);
        }
    }

    fn find_non_integral(&self, num_tests: usize) -> Option<u32> {
        let mut most_frac = None;
        for t in 0..num_tests {
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

    fn solve_lp(&self, node: NodeId, num_tests: usize) -> Option<SearchHead> {
        match self.lp_solver.solve() {
            INFEASIBLE => None,
            OPTIMAL => {
                let lb = NotNan::new(self.lp_solver.objective_value()).unwrap();
                let status = match self.find_non_integral(num_tests) {
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
}

fn collect_fixings(nodes: &[Node], mut node: NodeId) -> Vec<(u32, u8)> {
    let mut fixings = Vec::new();
    while node != NO_PARENT {
        if let Some(f) = nodes[node].fixing {
            fixings.push((f.var, f.val));
        }
        node = nodes[node].parent;
    }
    fixings
}

fn worker_loop(shared: &SharedState, worker: &WorkerSolver, num_tests: usize) {
    loop {
        let head = {
            let mut heads = shared.heads.lock().unwrap();
            loop {
                if shared.done.load(AtomicOrd::SeqCst) {
                    return;
                }
                if let Some(h) = heads.pop() {
                    shared.active.fetch_add(1, AtomicOrd::SeqCst);
                    break h;
                }
                // heap empty... are we done?
                if shared.active.load(AtomicOrd::SeqCst) == 0 {
                    shared.done.store(true, AtomicOrd::SeqCst);
                    shared.heads_cv.notify_all();
                    return;
                }
                // wait for more work
                heads = shared.heads_cv.wait(heads).unwrap();
            }
        };

        let cur_opt = *shared.ip_opt.lock().unwrap();
        if head.lb >= cur_opt {
            shared.active.fetch_sub(1, AtomicOrd::SeqCst);
            shared.heads_cv.notify_all();
            continue;
        }

        match head.status {
            NodeStatus::Integral => {
                let mut opt = shared.ip_opt.lock().unwrap();
                if head.lb < *opt {
                    *opt = head.lb;
                }
            }
            NodeStatus::Fractional { branch_var } => {
                // collect parent fixings
                let nodes = shared.nodes.read().unwrap();
                let parent_fixings = collect_fixings(&nodes, head.node);

                // insert child nodes
                let (l_id, r_id) = {
                    let mut nodes = shared.nodes.write().unwrap();
                    let l = nodes.len();
                    nodes.push(Node {
                        parent: head.node,
                        fixing: Some(Fixing {
                            var: branch_var,
                            val: 0,
                        }),
                    });
                    let r = nodes.len();
                    nodes.push(Node {
                        parent: head.node,
                        fixing: Some(Fixing {
                            var: branch_var,
                            val: 1,
                        }),
                    });
                    (l, r)
                };

                // solve left
                let mut fixings = parent_fixings.clone();
                fixings.push((branch_var, 0));
                worker.reset_and_apply(&fixings);
                worker.lp_solver.restore_basis(&head.basis);
                if let Some(l_head) = worker.solve_lp(l_id, num_tests) {
                    let mut heads = shared.heads.lock().unwrap();
                    heads.push(l_head);
                    shared.heads_cv.notify_one();
                }

                // solve right
                *fixings.last_mut().unwrap() = (branch_var, 1);
                worker.reset_and_apply(&fixings);
                worker.lp_solver.restore_basis(&head.basis);
                if let Some(r_head) = worker.solve_lp(r_id, num_tests) {
                    let mut heads = shared.heads.lock().unwrap();
                    heads.push(r_head);
                    shared.heads_cv.notify_one();
                }
            }
        }

        shared.active.fetch_sub(1, AtomicOrd::SeqCst);
        shared.heads_cv.notify_all();
    }
}

impl BNCSolver {
    pub fn new(instance: IPInstance) -> Self {
        BNCSolver { instance }
    }

    pub fn solve(&mut self) -> usize {
        let num_tests = self.instance.num_tests;
        let num_threads = 4;

        let root_worker = WorkerSolver::new(&self.instance);
        root_worker.lp_solver.solve();

        let Some(branch_var) = root_worker.find_non_integral(num_tests) else {
            return root_worker.lp_solver.objective_value() as usize;
        };

        let root_lb = NotNan::new(root_worker.lp_solver.objective_value()).unwrap();
        let root_basis = root_worker.lp_solver.save_basis();
        drop(root_worker);

        let shared = SharedState {
            nodes: RwLock::new(vec![Node {
                parent: NO_PARENT,
                fixing: None,
            }]),
            heads: Mutex::new(BinaryHeap::new()),
            heads_cv: Condvar::new(),
            ip_opt: Mutex::new(NotNan::new(f64::INFINITY).unwrap()),
            active: AtomicUsize::new(0),
            done: AtomicBool::new(false),
        };

        shared.heads.lock().unwrap().push(SearchHead {
            node: 0,
            lb: root_lb,
            status: NodeStatus::Fractional { branch_var },
            basis: root_basis,
        });

        let instance = &self.instance;
        thread::scope(|s| {
            for _ in 0..num_threads {
                s.spawn(|| {
                    let worker = WorkerSolver::new(instance);
                    worker_loop(&shared, &worker, num_tests);
                });
            }
        });

        let ip_opt = *shared.ip_opt.lock().unwrap();
        *ip_opt as usize
    }
}
