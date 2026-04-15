#include "ipinstance.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ortools/linear_solver/linear_solver.h"

namespace {

using operations_research::MPConstraint;
using operations_research::MPSolver;
using operations_research::MPVariable;

bool NearlyInteger(double value, double eps) {
  return std::abs(value - std::round(value)) <= eps;
}

}  // namespace

IPInstance::IPInstance(const std::string& filename)
    : best_objective_(std::numeric_limits<double>::infinity()) {
  LoadFromFile(filename);
  PreprocessInstance();
}

SolveResult IPInstance::Solve() {
  if (is_infeasible_) {
    return {};
  }

  const GreedyResult greedy = GreedyInitialSolution();
  if (!greedy.feasible) {
    return {};
  }

  best_solution_ = greedy.solution;
  best_objective_ = greedy.objective_value;
  has_best_solution_ = true;

  const std::vector<bool> fixed_zero(num_tests_, false);
  const std::vector<bool> fixed_one(num_tests_, false);
  BranchAndBound(fixed_zero, fixed_one);

  SolveResult result;
  result.has_solution = has_best_solution_;
  result.is_optimal = has_best_solution_;
  result.objective_value = best_objective_;
  result.solution = best_solution_;
  return result;
}

void IPInstance::LoadFromFile(const std::string& filename) {
  std::ifstream input(filename);
  if (!input) {
    throw std::runtime_error("Error reading instance file. File format may be incorrect.");
  }

  input >> num_tests_;
  input >> num_diseases_;

  if (!input || num_tests_ < 0 || num_diseases_ < 0) {
    throw std::runtime_error("Error reading instance file. File format may be incorrect.");
  }

  cost_of_test_.assign(num_tests_, 0.0);
  for (int i = 0; i < num_tests_; ++i) {
    input >> cost_of_test_[i];
  }

  a_.assign(num_tests_, std::vector<int>(num_diseases_, 0));
  for (int i = 0; i < num_tests_; ++i) {
    for (int j = 0; j < num_diseases_; ++j) {
      input >> a_[i][j];
    }
  }

  if (!input) {
    throw std::runtime_error("Error reading instance file. File format may be incorrect.");
  }
}

void IPInstance::PreprocessInstance() {
  pair_constraints_ = GetDistinguishingTestSets();
  for (const PairConstraint& pair : pair_constraints_) {
    if (pair.distinguishing_tests.empty()) {
      is_infeasible_ = true;
      return;
    }
  }

  is_infeasible_ = false;
  test_to_pairs_ = BuildTestPairCoverage();

  std::vector<int> keep;
  for (int i = 0; i < num_tests_; ++i) {
    if (!test_to_pairs_[i].empty()) {
      keep.push_back(i);
    }
  }

  std::vector<bool> dominated(num_tests_, false);
  for (int i : keep) {
    if (dominated[i]) {
      continue;
    }
    for (int j : keep) {
      if (i == j || dominated[j]) {
        continue;
      }
      if (cost_of_test_[i] <= cost_of_test_[j] &&
          IsSubset(test_to_pairs_[j], test_to_pairs_[i])) {
        dominated[j] = true;
      }
    }
  }

  std::vector<int> keep_list;
  keep_list.reserve(keep.size());
  for (int idx : keep) {
    if (!dominated[idx]) {
      keep_list.push_back(idx);
    }
  }

  if (static_cast<int>(keep_list.size()) == num_tests_) {
    return;
  }

  std::vector<double> reduced_costs;
  std::vector<std::vector<int>> reduced_a;
  reduced_costs.reserve(keep_list.size());
  reduced_a.reserve(keep_list.size());
  for (int idx : keep_list) {
    reduced_costs.push_back(cost_of_test_[idx]);
    reduced_a.push_back(a_[idx]);
  }

  cost_of_test_ = std::move(reduced_costs);
  a_ = std::move(reduced_a);
  num_tests_ = static_cast<int>(keep_list.size());

  pair_constraints_ = GetDistinguishingTestSets();
  test_to_pairs_ = BuildTestPairCoverage();
}

std::vector<PairConstraint> IPInstance::GetDistinguishingTestSets() const {
  std::vector<PairConstraint> constraints;
  for (int j = 0; j < num_diseases_; ++j) {
    for (int k = j + 1; k < num_diseases_; ++k) {
      PairConstraint pair{j, k, {}};
      for (int i = 0; i < num_tests_; ++i) {
        if (a_[i][j] != a_[i][k]) {
          pair.distinguishing_tests.push_back(i);
        }
      }
      constraints.push_back(std::move(pair));
    }
  }
  return constraints;
}

std::vector<std::vector<int>> IPInstance::BuildTestPairCoverage() const {
  std::vector<std::vector<int>> test_to_pairs(num_tests_);
  for (std::size_t pair_id = 0; pair_id < pair_constraints_.size(); ++pair_id) {
    for (int test_id : pair_constraints_[pair_id].distinguishing_tests) {
      test_to_pairs[test_id].push_back(static_cast<int>(pair_id));
    }
  }
  return test_to_pairs;
}

LPResult IPInstance::SolveLPRelaxation(const std::vector<bool>& fixed_zero,
                                       const std::vector<bool>& fixed_one) const {
  for (int i = 0; i < num_tests_; ++i) {
    if (fixed_zero[i] && fixed_one[i]) {
      return {};
    }
  }

  MPSolver solver("healthcare_lp", MPSolver::GLOP_LINEAR_PROGRAMMING);

  std::vector<MPVariable*> x;
  x.reserve(num_tests_);
  for (int i = 0; i < num_tests_; ++i) {
    double lb = 0.0;
    double ub = 1.0;
    if (fixed_zero[i]) {
      ub = 0.0;
    } else if (fixed_one[i]) {
      lb = 1.0;
    }
    x.push_back(solver.MakeNumVar(lb, ub, "x_" + std::to_string(i)));
  }

  auto* objective = solver.MutableObjective();
  for (int i = 0; i < num_tests_; ++i) {
    objective->SetCoefficient(x[i], cost_of_test_[i]);
  }
  objective->SetMinimization();

  for (const PairConstraint& pair : pair_constraints_) {
    if (pair.distinguishing_tests.empty()) {
      return {};
    }
    MPConstraint* constraint = solver.MakeRowConstraint(1.0, solver.infinity());
    for (int test_id : pair.distinguishing_tests) {
      constraint->SetCoefficient(x[test_id], 1.0);
    }
  }

  const MPSolver::ResultStatus status = solver.Solve();
  if (status != MPSolver::OPTIMAL) {
    return {};
  }

  LPResult result;
  result.feasible = true;
  result.objective_value = objective->Value();
  result.x_values.resize(num_tests_);
  result.is_integral = true;
  for (int i = 0; i < num_tests_; ++i) {
    result.x_values[i] = x[i]->solution_value();
    if (!NearlyInteger(result.x_values[i], kEps)) {
      result.is_integral = false;
    }
  }
  return result;
}

int IPInstance::ChooseBranchVariable(const std::vector<double>& x_values,
                                     const std::vector<bool>& fixed_zero,
                                     const std::vector<bool>& fixed_one) const {
  int best_var = -1;
  int best_coverage = -1;
  double best_distance_to_half = std::numeric_limits<double>::infinity();

  for (int i = 0; i < num_tests_; ++i) {
    if (fixed_zero[i] || fixed_one[i] || NearlyInteger(x_values[i], kEps)) {
      continue;
    }

    const int coverage = static_cast<int>(test_to_pairs_[i].size());
    const double distance_to_half = std::abs(x_values[i] - 0.5);
    if (coverage > best_coverage ||
        (coverage == best_coverage && distance_to_half < best_distance_to_half)) {
      best_var = i;
      best_coverage = coverage;
      best_distance_to_half = distance_to_half;
    }
  }

  return best_var;
}

std::vector<int> IPInstance::ExtractIntegerSolution(const std::vector<double>& x_values) const {
  std::vector<int> solution(x_values.size(), 0);
  for (std::size_t i = 0; i < x_values.size(); ++i) {
    solution[i] = static_cast<int>(std::llround(x_values[i]));
  }
  return solution;
}

std::vector<int> IPInstance::RemoveRedundantTests(const std::vector<int>& chosen_tests) const {
  std::vector<int> coverage_count(pair_constraints_.size(), 0);
  for (int test_id : chosen_tests) {
    for (int pair_id : test_to_pairs_[test_id]) {
      ++coverage_count[pair_id];
    }
  }

  std::vector<int> ordered = chosen_tests;
  std::sort(ordered.begin(), ordered.end(), [&](int lhs, int rhs) {
    return cost_of_test_[lhs] > cost_of_test_[rhs];
  });

  std::vector<bool> removed(num_tests_, false);
  for (int test_id : ordered) {
    bool can_remove = true;
    for (int pair_id : test_to_pairs_[test_id]) {
      if (coverage_count[pair_id] <= 1) {
        can_remove = false;
        break;
      }
    }
    if (!can_remove) {
      continue;
    }
    removed[test_id] = true;
    for (int pair_id : test_to_pairs_[test_id]) {
      --coverage_count[pair_id];
    }
  }

  std::vector<int> remaining;
  for (int test_id : chosen_tests) {
    if (!removed[test_id]) {
      remaining.push_back(test_id);
    }
  }
  return remaining;
}

IPInstance::GreedyResult IPInstance::BuildFeasibleSolution(
  const std::vector<int>& initial_chosen_tests,
  const std::vector<bool>& excluded_tests,
  const std::function<double(int, const std::vector<int>&)>& score_fn) const {
  GreedyResult result;

  std::vector<bool> chosen(num_tests_, false);
  std::vector<bool> uncovered(pair_constraints_.size(), true);
  int uncovered_count = static_cast<int>(pair_constraints_.size());

  for (int test_id : initial_chosen_tests) {
    if (test_id < 0 || test_id >= num_tests_ || excluded_tests[test_id]) {
      return result;
    }
    if (chosen[test_id]) {
      continue;
    }
    chosen[test_id] = true;
    for (int pair_id : test_to_pairs_[test_id]) {
      if (uncovered[pair_id]) {
        uncovered[pair_id] = false;
        --uncovered_count;
      }
    }
  }

  while (uncovered_count > 0) {
    int best_test = -1;
    double best_score = -1.0;
    std::vector<int> best_new_pairs;

    for (int i = 0; i < num_tests_; ++i) {
      if (chosen[i] || excluded_tests[i]) {
        continue;
      }

      std::vector<int> new_pairs;
      for (int pair_id : test_to_pairs_[i]) {
        if (uncovered[pair_id]) {
          new_pairs.push_back(pair_id);
        }
      }

      if (new_pairs.empty()) {
        continue;
      }

      const double score = score_fn(i, new_pairs);
      if (score > best_score) {
        best_test = i;
        best_score = score;
        best_new_pairs = std::move(new_pairs);
      }
    }

    if (best_test == -1) {
      return result;
    }

    chosen[best_test] = true;
    for (int pair_id : best_new_pairs) {
      if (uncovered[pair_id]) {
        uncovered[pair_id] = false;
        --uncovered_count;
      }
    }
  }

  std::vector<int> chosen_tests;
  for (int i = 0; i < num_tests_; ++i) {
    if (chosen[i]) {
      chosen_tests.push_back(i);
    }
  }

  chosen_tests = RemoveRedundantTests(chosen_tests);

  result.solution.assign(num_tests_, 0);
  result.objective_value = 0.0;
  for (int test_id : chosen_tests) {
    result.solution[test_id] = 1;
    result.objective_value += cost_of_test_[test_id];
  }
  result.feasible = true;

  return result;
}


IPInstance::GreedyResult IPInstance::GreedyInitialSolution() const {
  const std::vector<int> initial_chosen_tests;
  const std::vector<bool> excluded_tests(num_tests_, false);

  return BuildFeasibleSolution(
      initial_chosen_tests,
      excluded_tests,
      [this](int test_id, const std::vector<int>& new_pairs) {
        return static_cast<double>(new_pairs.size()) / cost_of_test_[test_id];
      });
}

IPInstance::GreedyResult IPInstance::HeuristicSolutionFromLP(
  const std::vector<double>& x_values,
  const std::vector<bool>& fixed_zero,
  const std::vector<bool>& fixed_one) const {
  std::vector<int> initial_chosen_tests;
  for (int i = 0; i < num_tests_; ++i) {
    if (fixed_one[i]) {
      initial_chosen_tests.push_back(i);
    }
  }

  return BuildFeasibleSolution(
      initial_chosen_tests,
      fixed_zero,
      [this, &x_values](int test_id, const std::vector<int>& new_pairs) {
        return x_values[test_id] +
              static_cast<double>(new_pairs.size()) / cost_of_test_[test_id];
      });
}


double IPInstance::ForcedOneCost(const std::vector<bool>& fixed_one) const {
  double total = 0.0;
  for (int i = 0; i < num_tests_; ++i) {
    if (fixed_one[i]) {
      total += cost_of_test_[i];
    }
  }
  return total;
}

void IPInstance::BranchAndBound(const std::vector<bool>& fixed_zero,
                                const std::vector<bool>& fixed_one) {
  if (ForcedOneCost(fixed_one) >= best_objective_ - kEps) {
    return;
  }

  const LPResult result = SolveLPRelaxation(fixed_zero, fixed_one);
  if (!result.feasible) {
    return;
  }

  if (result.objective_value >= best_objective_ - kEps) {
    return;
  }

  if (result.is_integral) {
    best_objective_ = result.objective_value;
    best_solution_ = ExtractIntegerSolution(result.x_values);
    has_best_solution_ = true;
    return;
  }

  const GreedyResult heuristic = HeuristicSolutionFromLP(
    result.x_values, fixed_zero, fixed_one);
    
  if (heuristic.feasible && heuristic.objective_value < best_objective_) {
    best_objective_ = heuristic.objective_value;
    best_solution_ = heuristic.solution;
    has_best_solution_ = true;
  }

  const int branch_var = ChooseBranchVariable(result.x_values, fixed_zero, fixed_one);
  if (branch_var == -1) {
    return;
  }

  std::vector<bool> left_zero = fixed_zero;
  left_zero[branch_var] = true;
  BranchAndBound(left_zero, fixed_one);

  std::vector<bool> right_one = fixed_one;
  right_one[branch_var] = true;
  BranchAndBound(fixed_zero, right_one);
}

bool IPInstance::IsSubset(const std::vector<int>& subset, const std::vector<int>& superset) {
  return std::includes(superset.begin(), superset.end(), subset.begin(), subset.end());
}
