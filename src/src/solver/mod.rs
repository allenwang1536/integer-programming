use std::{
    cmp::Ordering,
    collections::{BinaryHeap, HashSet},
    env, f64,
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering as AtomicOrd},
        Condvar, Mutex, RwLock,
    },
    thread,
};

struct FixingState {
    fixed: Vec<u64>,
    values: Vec<u64>,
}

impl FixingState {
    fn new(num_vars: usize) -> Self {
        let words = (num_vars + 63) / 64;
        FixingState {
            fixed: vec![0u64; words],
            values: vec![0u64; words],
        }
    }

    fn apply_transition(&mut self, new_fixings: &[(u32, u8)], worker: &WorkerSolver) {
        let mut new_fixed = vec![0u64; self.fixed.len()];
        let mut new_values = vec![0u64; self.fixed.len()];

        for &(var, val) in new_fixings {
            let w = var as usize / 64;
            let b = var as u64 % 64;
            new_fixed[w] |= 1 << b;
            if val == 1 {
                new_values[w] |= 1 << b;
            }
        }

        for w in 0..self.fixed.len() {
            let fix_changed = self.fixed[w] ^ new_fixed[w];
            let val_changed = self.values[w] ^ new_values[w];
            let changed = fix_changed | (new_fixed[w] & val_changed);

            let mut bits = changed;
            while bits != 0 {
                let bit = bits.trailing_zeros() as usize;
                let var = w * 64 + bit;
                let is_fixed = (new_fixed[w] >> bit) & 1;
                if is_fixed == 0 {
                    worker.test_var[var].set_bounds(0.0, 1.0);
                } else {
                    let val = ((new_values[w] >> bit) & 1) as f64;
                    worker.test_var[var].set_bounds(val, val);
                }
                bits &= bits - 1;
            }
        }

        self.fixed = new_fixed;
        self.values = new_values;
    }
}

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

#[derive(Clone)]
struct CutRow {
    vars: Vec<u32>,
    lb: f64,
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
    explored: AtomicUsize,
    done: AtomicBool,
}

struct WorkerSolver {
    lp_solver: OrSolverHandle,
    test_var: Vec<OrVarHandle>,
}

#[derive(Default)]
struct RootCutStats {
    initial_lb: Option<f64>,
    final_lb: f64,
    cuts_added: usize,
    lp_solves: usize,
}

impl WorkerSolver {
    fn new(instance: &IPInstance, cuts: &[CutRow]) -> Self {
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

        let worker = WorkerSolver {
            lp_solver,
            test_var,
        };
        for cut in cuts {
            worker.add_cut(cut);
        }

        worker
    }

    fn add_cut(&self, cut: &CutRow) {
        let c = self.lp_solver.make_constraint(cut.lb, f64::INFINITY);
        for &var in &cut.vars {
            self.lp_solver
                .constraint_set_coefficient(&c, &self.test_var[var as usize], 1.0);
        }
    }

    fn solution_values(&self, num_tests: usize) -> Vec<f64> {
        let mut values = Vec::with_capacity(num_tests);
        for t in 0..num_tests {
            values.push(self.test_var[t].solution_value());
        }
        values
    }

    fn run_root_cuts(&self, instance: &IPInstance) -> (Vec<CutRow>, RootCutStats) {
        const MAX_CUT_PASSES: usize = 3;
        const MAX_TRIPLET_CUTS_PER_PASS: usize = 64;

        let mut cuts = Vec::new();
        let mut stats = RootCutStats::default();
        let mut seen = HashSet::<Vec<u32>>::new();

        for _ in 0..MAX_CUT_PASSES {
            match self.lp_solver.solve() {
                INFEASIBLE => {
                    stats.final_lb = f64::INFINITY;
                    return (cuts, stats);
                }
                OPTIMAL => {}
                _ => unreachable!(),
            }

            stats.lp_solves += 1;
            let lb = self.lp_solver.objective_value();
            stats.initial_lb.get_or_insert(lb);

            let solution = self.solution_values(instance.num_tests);
            let mut new_cuts = Vec::new();

            if let Some(cut) = separate_global_counting_cut(instance, &solution) {
                push_new_cut(&mut new_cuts, &mut seen, cut);
            }

            for cut in separate_triplet_cuts(instance, &solution, MAX_TRIPLET_CUTS_PER_PASS) {
                push_new_cut(&mut new_cuts, &mut seen, cut);
            }

            if new_cuts.is_empty() {
                stats.final_lb = lb;
                return (cuts, stats);
            }

            for cut in &new_cuts {
                self.add_cut(cut);
            }
            stats.cuts_added += new_cuts.len();
            cuts.extend(new_cuts);
        }

        match self.lp_solver.solve() {
            INFEASIBLE => stats.final_lb = f64::INFINITY,
            OPTIMAL => {
                stats.lp_solves += 1;
                stats.final_lb = self.lp_solver.objective_value();
            }
            _ => unreachable!(),
        }

        (cuts, stats)
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

fn stats_enabled() -> bool {
    env::var("IP_STATS").ok().as_deref() == Some("1")
}

fn counting_rhs(size: usize) -> f64 {
    size.next_power_of_two().trailing_zeros() as f64
}

fn push_new_cut(new_cuts: &mut Vec<CutRow>, seen: &mut HashSet<Vec<u32>>, cut: CutRow) {
    if seen.insert(cut.vars.clone()) {
        new_cuts.push(cut);
    }
}

fn separate_global_counting_cut(instance: &IPInstance, solution: &[f64]) -> Option<CutRow> {
    let mut vars = Vec::new();
    let mut lhs = 0.0;

    for t in 0..instance.num_tests {
        let first = instance.at(t, 0);
        let splits_all = (1..instance.num_diseases).any(|d| instance.at(t, d) != first);
        if splits_all {
            vars.push(t as u32);
            lhs += solution[t];
        }
    }

    let rhs = counting_rhs(instance.num_diseases);
    (lhs + 1e-9 < rhs).then_some(CutRow { vars, lb: rhs })
}

fn separate_triplet_cuts(instance: &IPInstance, solution: &[f64], max_cuts: usize) -> Vec<CutRow> {
    let mut violated = Vec::<(f64, Vec<u32>)>::new();

    for a in 0..instance.num_diseases {
        for b in a + 1..instance.num_diseases {
            for c in b + 1..instance.num_diseases {
                let mut vars = Vec::new();
                let mut lhs = 0.0;
                for t in 0..instance.num_tests {
                    let v0 = instance.at(t, a);
                    if instance.at(t, b) != v0 || instance.at(t, c) != v0 {
                        vars.push(t as u32);
                        lhs += solution[t];
                    }
                }

                if lhs + 1e-9 >= 2.0 || vars.is_empty() {
                    continue;
                }

                violated.push((lhs, vars));
            }
        }
    }

    violated.sort_by(|(lhs_a, _), (lhs_b, _)| lhs_a.total_cmp(lhs_b));
    violated
        .into_iter()
        .take(max_cuts)
        .map(|(_, vars)| CutRow { vars, lb: 2.0 })
        .collect()
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
    let mut fixing_state = FixingState::new(num_tests);
    loop {
        let head = {
            let mut heads = shared.heads.lock().unwrap();
            loop {
                if shared.done.load(AtomicOrd::SeqCst) {
                    return;
                }
                if let Some(h) = heads.pop() {
                    shared.active.fetch_add(1, AtomicOrd::SeqCst);
                    shared.explored.fetch_add(1, AtomicOrd::SeqCst);
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
                let parent_fixings = {
                    let nodes = shared.nodes.read().unwrap();
                    collect_fixings(&nodes, head.node)
                };

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
                fixing_state.apply_transition(&fixings, worker);
                worker.lp_solver.restore_basis(&head.basis);
                if let Some(l_head) = worker.solve_lp(l_id, num_tests) {
                    let mut heads = shared.heads.lock().unwrap();
                    heads.push(l_head);
                    shared.heads_cv.notify_one();
                }

                // solve right
                *fixings.last_mut().unwrap() = (branch_var, 1);
                fixing_state.apply_transition(&fixings, worker);
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
        let num_threads = 8;

        let root_worker = WorkerSolver::new(&self.instance, &[]);
        let (root_cuts, root_cut_stats) = root_worker.run_root_cuts(&self.instance);

        let Some(branch_var) = root_worker.find_non_integral(num_tests) else {
            if stats_enabled() {
                eprintln!(
                    "cuts=global+triplets root_lb_before={:.6} root_lb_after={:.6} root_cuts={} root_lp_solves={} nodes=0",
                    root_cut_stats.initial_lb.unwrap_or(0.0),
                    root_cut_stats.final_lb,
                    root_cut_stats.cuts_added,
                    root_cut_stats.lp_solves
                );
            }
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
            explored: AtomicUsize::new(0),
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
                    let worker = WorkerSolver::new(instance, &root_cuts);
                    worker_loop(&shared, &worker, num_tests);
                });
            }
        });

        let ip_opt = *shared.ip_opt.lock().unwrap();
        if stats_enabled() {
            eprintln!(
                "cuts=global+triplets root_lb_before={:.6} root_lb_after={:.6} root_cuts={} root_lp_solves={} nodes={}",
                root_cut_stats.initial_lb.unwrap_or(0.0),
                root_cut_stats.final_lb,
                root_cut_stats.cuts_added,
                root_cut_stats.lp_solves,
                shared.explored.load(AtomicOrd::SeqCst)
            );
        }
        *ip_opt as usize
    }
}
