#include <cstdint>
#include <cstring>
#include <vector>

#include "ortools/linear_solver/linear_solver.h"

namespace or_tools = operations_research;

struct LpBasis {
  std::vector<double> var_values;
};

extern "C" {

typedef void* OrSolverPtr;
typedef void* OrVarPtr;
typedef void* OrConstraintPtr;
typedef void* OrBasisPtr;

OrSolverPtr or_new_mpsolver_scip() {
  auto* solver = new or_tools::MPSolver(
      "SCIP",
      or_tools::MPSolver::SCIP_MIXED_INTEGER_PROGRAMMING);
  return reinterpret_cast<OrSolverPtr>(solver);
}

extern "C" OrVarPtr or_mpsolver_make_num_var(OrSolverPtr p, double lb, double ub, const char* name) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  auto* var = solver->MakeNumVar(lb, ub, name ? std::string(name) : "");
  return reinterpret_cast<OrVarPtr>(var);
}

extern "C" OrVarPtr or_mpsolver_make_int_var(OrSolverPtr p, double lb, double ub, const char* name) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  auto* var = solver->MakeIntVar(lb, ub, name ? std::string(name) : "");
  return reinterpret_cast<OrVarPtr>(var);
}

extern "C" OrConstraintPtr or_mpsolver_make_constraint(OrSolverPtr p, double lb, double ub) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  auto* ct = solver->MakeRowConstraint(lb, ub);
  return reinterpret_cast<OrConstraintPtr>(ct);
}

extern "C" void or_constraint_set_coefficient(OrConstraintPtr c, OrVarPtr v, double coeff) {
  auto* ct = reinterpret_cast<or_tools::MPConstraint*>(c);
  auto* var = reinterpret_cast<or_tools::MPVariable*>(v);
  ct->SetCoefficient(var, coeff);
}

extern "C" void or_objective_set_coefficient(OrSolverPtr p, OrVarPtr v, double coeff) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  auto* var = reinterpret_cast<or_tools::MPVariable*>(v);
  solver->MutableObjective()->SetCoefficient(var, coeff);
}

extern "C" void or_objective_set_maximize(OrSolverPtr p) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  solver->MutableObjective()->SetMaximization();
}

extern "C" void or_objective_set_minimize(OrSolverPtr p) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  solver->MutableObjective()->SetMinimization();
}

extern "C" void or_mpsolver_set_time_limit_ms(OrSolverPtr p, long long ms) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  solver->set_time_limit(ms);
}

extern "C" int or_mpsolver_solve(OrSolverPtr p) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  return static_cast<int>(solver->Solve());
}

extern "C" OrBasisPtr or_mpsolver_save_basis(OrSolverPtr p) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  auto* basis = new LpBasis();
  const auto& vars = solver->variables();
  basis->var_values.resize(vars.size());
  for (size_t i = 0; i < vars.size(); ++i) {
    basis->var_values[i] = vars[i]->solution_value();
  }
  return reinterpret_cast<OrBasisPtr>(basis);
}

extern "C" void or_mpsolver_restore_basis(OrSolverPtr p, OrBasisPtr b) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  auto* basis = reinterpret_cast<LpBasis*>(b);
  const auto& vars = solver->variables();
  std::vector<std::pair<const or_tools::MPVariable*, double>> hint;
  hint.reserve(vars.size());
  for (size_t i = 0; i < vars.size(); ++i) {
    hint.emplace_back(vars[i], basis->var_values[i]);
  }
  solver->SetHint(std::move(hint));
}

extern "C" void or_delete_basis(OrBasisPtr b) {
  delete reinterpret_cast<LpBasis*>(b);
}

extern "C" double or_mpsolver_objective_value(OrSolverPtr p) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  return solver->Objective().Value();
}

extern "C" double or_var_solution_value(OrVarPtr v) {
  auto* var = reinterpret_cast<or_tools::MPVariable*>(v);
  return var->solution_value();
}

extern "C" void or_var_set_bounds(OrVarPtr v, double lb, double ub) {
  auto* var = reinterpret_cast<or_tools::MPVariable*>(v);
  var->SetBounds(lb, ub);
}

void or_delete_mpsolver(OrSolverPtr p) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  delete solver;
}

int32_t or_mpsolver_version(OrSolverPtr p, char* buf, int32_t buf_len) {
  auto* solver = reinterpret_cast<or_tools::MPSolver*>(p);
  const std::string ver = solver->SolverVersion();
  const int32_t need = static_cast<int32_t>(ver.size());
  if (buf && buf_len > 0) {
    const int32_t copy_len = need < buf_len ? need : buf_len - 1;
    if (copy_len > 0) {
      std::memcpy(buf, ver.data(), static_cast<size_t>(copy_len));
    }
    buf[copy_len] = '\0';
  }
  return need;
}

}
