#ifndef IPINSTANCE_H
#define IPINSTANCE_H

#include <cstddef>
#include <string>
#include <vector>
#include <functional>


struct PairConstraint {
  int disease_a;
  int disease_b;
  std::vector<int> distinguishing_tests;
};

struct LPResult {
  bool feasible = false;
  double objective_value = 0.0;
  std::vector<double> x_values;
  bool is_integral = false;
};

struct SolveResult {
  bool has_solution = false;
  bool is_optimal = false;
  double objective_value = 0.0;
  std::vector<int> solution;
};

class IPInstance {
 public:
  explicit IPInstance(const std::string& filename);

  SolveResult Solve();

 private:
  struct GreedyResult {
    bool feasible = false;
    double objective_value = 0.0;
    std::vector<int> solution;
  };

  void LoadFromFile(const std::string& filename);
  void PreprocessInstance();

  std::vector<PairConstraint> GetDistinguishingTestSets() const;
  std::vector<std::vector<int>> BuildTestPairCoverage() const;

  LPResult SolveLPRelaxation(const std::vector<bool>& fixed_zero,
                             const std::vector<bool>& fixed_one) const;
  int ChooseBranchVariable(const std::vector<double>& x_values,
                           const std::vector<bool>& fixed_zero,
                           const std::vector<bool>& fixed_one) const;
  std::vector<int> ExtractIntegerSolution(const std::vector<double>& x_values) const;

  GreedyResult GreedyInitialSolution() const;
  GreedyResult BuildFeasibleSolution(
      const std::vector<int>& initial_chosen_tests,
      const std::vector<bool>& excluded_tests,
      const std::function<double(int, const std::vector<int>&)>& score_fn) const;

  GreedyResult HeuristicSolutionFromLP(const std::vector<double>& x_values,
                                       const std::vector<bool>& fixed_zero,
                                       const std::vector<bool>& fixed_one) const;
  std::vector<int> RemoveRedundantTests(const std::vector<int>& chosen_tests) const;

  bool PropagateForcedChoices(std::vector<bool>& fixed_zero,
                              std::vector<bool>& fixed_one) const;
  double ForcedOneCost(const std::vector<bool>& fixed_one) const;
  void BranchAndBound(const std::vector<bool>& fixed_zero,
                      const std::vector<bool>& fixed_one);

  static bool IsSubset(const std::vector<int>& subset, const std::vector<int>& superset);

  int num_tests_ = 0;
  int num_diseases_ = 0;
  std::vector<double> cost_of_test_;
  std::vector<std::vector<int>> a_;

  std::vector<PairConstraint> pair_constraints_;
  std::vector<std::vector<int>> test_to_pairs_;

  bool is_infeasible_ = false;
  bool has_best_solution_ = false;
  double best_objective_ = 0.0;
  std::vector<int> best_solution_;

  static constexpr double kEps = 1e-6;
};

#endif
