use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_longlong, c_void};

extern "C" {
    fn or_new_mpsolver_scip() -> *mut c_void;
    fn or_delete_mpsolver(ptr: *mut c_void);
    fn or_mpsolver_version(ptr: *mut c_void, buf: *mut c_char, buf_len: c_int) -> c_int;
    fn or_mpsolver_make_num_var(
        ptr: *mut c_void,
        lb: f64,
        ub: f64,
        name: *const c_char,
    ) -> *mut c_void;
    fn or_mpsolver_make_int_var(
        ptr: *mut c_void,
        lb: f64,
        ub: f64,
        name: *const c_char,
    ) -> *mut c_void;
    fn or_mpsolver_make_constraint(ptr: *mut c_void, lb: f64, ub: f64) -> *mut c_void;
    fn or_constraint_set_coefficient(ct: *mut c_void, var: *mut c_void, coeff: f64);
    fn or_objective_set_coefficient(ptr: *mut c_void, var: *mut c_void, coeff: f64);
    fn or_objective_set_maximize(ptr: *mut c_void);
    fn or_objective_set_minimize(ptr: *mut c_void);
    fn or_mpsolver_set_time_limit_ms(ptr: *mut c_void, ms: c_longlong);
    fn or_mpsolver_solve(ptr: *mut c_void) -> c_int;
    fn or_mpsolver_objective_value(ptr: *mut c_void) -> f64;
    fn or_var_solution_value(var: *mut c_void) -> f64;
    fn or_var_set_bounds(ptr: *mut c_void, lb: f64, ub: f64);
}

#[repr(transparent)]
#[derive(Clone)]
pub struct OrSolverHandle(pub *mut c_void);

#[repr(transparent)]
#[derive(Clone)]
pub struct OrVarHandle(pub *mut c_void);

#[repr(transparent)]
#[derive(Clone)]
pub struct OrConstraintHandle(pub *mut c_void);

impl OrSolverHandle {
    pub fn new_scip() -> Self {
        let p = unsafe { or_new_mpsolver_scip() };
        Self(p)
    }
    pub fn version(&self) -> String {
        let need = unsafe { or_mpsolver_version(self.0, std::ptr::null_mut(), 0) } as usize;
        let mut buf = vec![0u8; need + 1];
        unsafe {
            or_mpsolver_version(self.0, buf.as_mut_ptr() as *mut c_char, buf.len() as c_int);
        }
        let cstr = unsafe { CStr::from_ptr(buf.as_ptr() as *const c_char) };
        cstr.to_string_lossy().into_owned()
    }

    pub fn make_num_var(&self, lb: f64, ub: f64, name: &str) -> OrVarHandle {
        let cname = CString::new(name).unwrap_or_default();
        let p = unsafe { or_mpsolver_make_num_var(self.0, lb, ub, cname.as_ptr()) };
        OrVarHandle(p)
    }
    pub fn make_int_var(&self, lb: f64, ub: f64, name: &str) -> OrVarHandle {
        let cname = CString::new(name).unwrap_or_default();
        let p = unsafe { or_mpsolver_make_int_var(self.0, lb, ub, cname.as_ptr()) };
        OrVarHandle(p)
    }
    pub fn make_constraint(&self, lb: f64, ub: f64) -> OrConstraintHandle {
        let p = unsafe { or_mpsolver_make_constraint(self.0, lb, ub) };
        OrConstraintHandle(p)
    }
    pub fn constraint_set_coefficient(
        &self,
        ct: &OrConstraintHandle,
        var: &OrVarHandle,
        coeff: f64,
    ) {
        unsafe { or_constraint_set_coefficient(ct.0, var.0, coeff) }
    }
    pub fn objective_set_coefficient(&self, var: &OrVarHandle, coeff: f64) {
        unsafe { or_objective_set_coefficient(self.0, var.0, coeff) }
    }
    pub fn objective_set_maximize(&self) {
        unsafe { or_objective_set_maximize(self.0) }
    }
    pub fn objective_set_minimize(&self) {
        unsafe { or_objective_set_minimize(self.0) }
    }
    pub fn set_time_limit_ms(&self, ms: i64) {
        unsafe { or_mpsolver_set_time_limit_ms(self.0, ms as c_longlong) }
    }
    pub fn solve(&self) -> i32 {
        unsafe { or_mpsolver_solve(self.0) as i32 }
    }
    pub fn objective_value(&self) -> f64 {
        unsafe { or_mpsolver_objective_value(self.0) }
    }
}

impl OrVarHandle {
    pub fn solution_value(&self) -> f64 {
        unsafe { or_var_solution_value(self.0) }
    }

    pub fn set_bounds(&self, lb: f64, ub: f64) {
        unsafe { or_var_set_bounds(self.0, lb, ub) }
    }
}

impl Drop for OrSolverHandle {
    fn drop(&mut self) {
        unsafe { or_delete_mpsolver(self.0) }
    }
}

pub use OrConstraintHandle as ConstraintHandle;
pub use OrSolverHandle as SolverHandle;
pub use OrVarHandle as VarHandle;

// Mirror MPSolver::ResultStatus for convenience.
pub mod status {
    pub const OPTIMAL: i32 = 0;
    pub const FEASIBLE: i32 = 1;
    pub const INFEASIBLE: i32 = 2;
    pub const UNBOUNDED: i32 = 3;
    pub const ABNORMAL: i32 = 4;
    pub const MODEL_INVALID: i32 = 5;
    pub const NOT_SOLVED: i32 = 6;
}
