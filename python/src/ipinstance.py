from itertools import count

import heapq
import numpy as np
from ortools.linear_solver import pywraplp

#  * File Format
#  * #Tests (i.e., n)
#  * #Diseases (i.e., m)
#  * Cost_1 Cost_2 . . . Cost_n
#  * A(1,1) A(1,2) . . . A(1, m)
#  * A(2,1) A(2,2) . . . A(2, m)
#  * . . . . . . . . . . . . . .
#  * A(n,1) A(n,2) . . . A(n, m)

class IPInstance:
    numTests: int  # number of tests
    numDiseases: int  # number of diseases
    costOfTest: np.ndarray  # [numTests] the cost of each test
    A: np.ndarray  # [numTests][numDiseases] 0/1 matrix if test is positive for disease

    def __init__(self, filename: str) -> None:
        self.load_from_file(filename)
        self.solution = None
        self.objective_value = None

        self.best_objective = float("inf")
        self.best_solution = None
        self.is_infeasible = False

        self.preprocess_instance()

    def preprocess_instance(self):
        """
        1. detect infeasibility
        2. remove useless tests
        3. remove dominated tests
        """

        self.pair_constraints = self.get_distinguishing_test_sets()
        for _, _, distinguishing_tests in self.pair_constraints:
            if not distinguishing_tests:
                self.is_infeasible = True
                return

        self.is_infeasible = False
        self.test_to_pairs = self.build_test_pair_coverage()

        keep = {i for i in range(self.numTests) if len(self.test_to_pairs[i]) > 0}

        # keep = sorted(keep)
        dominated = set()

        for i in keep:
            if i in dominated:
                continue
                
            for j in keep:
                if i == j or j in dominated:
                    continue
                    
                if (self.test_to_pairs[j].issubset(self.test_to_pairs[i]) and 
                    self.costOfTest[i] <= self.costOfTest[j]):
                    dominated.add(j)
        
        keep -= dominated


        # if nothing changed just return
        if len(keep) == self.numTests:
            return

        keep_list = sorted(keep) # sort for determinism

        self.costOfTest = self.costOfTest[keep_list]
        self.A = self.A[keep_list, :]
        self.numTests = len(keep_list)

        # rebuild now since our input structures are now different
        self.pair_constraints = self.get_distinguishing_test_sets()
        self.test_to_pairs = self.build_test_pair_coverage()


    def get_distinguishing_test_sets(self):
        pair_constraints = []

        for j in range(self.numDiseases):
            for k in range(j + 1, self.numDiseases):
                distinguishing_tests = [
                    i for i in range(self.numTests)
                    if self.A[i, j] != self.A[i, k]
                ]
                pair_constraints.append((j, k, distinguishing_tests))

        return pair_constraints

    def solve_lp_relaxation(self, fixed_zero=None, fixed_one=None):
        # default arguments are shared across all function calls, so set to None and init

        if fixed_zero is None:
            fixed_zero = set()
        if fixed_one is None:
            fixed_one = set()

        # inconsistent branching decision
        if fixed_zero & fixed_one:
            return {
                "feasible": False,
                "objective_value": None,
                "x_values": None, 
                "is_integral": False,
            }
        
        solver = pywraplp.Solver.CreateSolver("SCIP")
        
        # decision variables: whether we should include / exclude tests
        x = []
        for i in range(self.numTests):
            if i in fixed_zero:
                var = solver.NumVar(0.0, 0.0, f"x_{i}")
            elif i in fixed_one:
                var = solver.NumVar(1.0, 1.0, f"x_{i}")
            else:
                var = solver.NumVar(0.0, 1.0, f"x_{i}")
            x.append(var)
        
        # objective: minimize total cost of tests
        objective = solver.Objective()
        for i in range(self.numTests):
            objective.SetCoefficient(x[i], float(self.costOfTest[i]))
        objective.SetMinimization()

        # constraints: each pair of diseases must be distinguishable
        for j, k, distinguishing_tests in self.pair_constraints:

            # infeasible since all tests is not a solution
            if not distinguishing_tests:
                return {
                    "feasible": False,
                    "objective_value": None,
                    "x_values": None,
                    "is_integral": False,
                }

            constraint = solver.RowConstraint(1.0, solver.infinity(), f"pair_{j}_{k}")
            for i in distinguishing_tests:
                constraint.SetCoefficient(x[i], 1.0)

        status = solver.Solve()
        if status != pywraplp.Solver.OPTIMAL:
            return {
                "feasible": False,
                "objective_value": None,
                "x_values": None,
                "is_integral": False,
            }
        
        x_values = [x[i].solution_value() for i in range(self.numTests)]
        eps = 1e-6
        is_integral = all(abs(v - round(v)) <= eps for v in x_values)

        return {
            "feasible": True,
            "objective_value": solver.Objective().Value(),
            "x_values": x_values,
            "is_integral": is_integral,
        }
    
    def choose_branch_variable(self, x_values, fixed_zero, fixed_one, eps=1e-6):
        best_var = None
        best_coverage = -1
        best_distance_to_half = float('inf')
        for i in range(self.numTests):
            if i in fixed_zero or i in fixed_one:
                continue
        
            if abs(x_values[i] - round(x_values[i])) <= eps:
                continue

            coverage = len(self.test_to_pairs[i])
            distance_to_half = abs(x_values[i] - 0.5)

            if (coverage > best_coverage or 
                (coverage == best_coverage and distance_to_half < best_distance_to_half)
                ):
                best_var = i
                best_coverage = coverage
                best_distance_to_half = distance_to_half
        
        return best_var


    def extract_integer_solution(self, x_values):
        return [round(v) for v in x_values]

    def build_test_pair_coverage(self):
        test_to_pairs = [set() for _ in range(self.numTests)]

        for pair_id, (_, _, distinguishing_tests) in enumerate(self.pair_constraints):
            for test_id in distinguishing_tests:
                test_to_pairs[test_id].add(pair_id)

        return test_to_pairs

    def remove_redundant_tests(self, chosen_tests):
        chosen_tests = set(chosen_tests)

        for test_id in sorted(chosen_tests, key=lambda i: self.costOfTest[i], reverse=True):
            trial_tests = chosen_tests - {test_id}

            covered_pairs = set()
            for i in trial_tests:
                covered_pairs |= self.test_to_pairs[i]

            if len(covered_pairs) == len(self.pair_constraints):
                chosen_tests.remove(test_id)

        return chosen_tests

    def build_feasible_solution(self, initial_chosen_tests, excluded_tests, score_fn):
        chosen_tests = set(initial_chosen_tests)
        excluded_tests = set(excluded_tests)

        uncovered_pairs = set(range(len(self.pair_constraints)))

        # remove pairs already covered by chosen tests
        for test_id in chosen_tests:
            uncovered_pairs -= self.test_to_pairs[test_id]

        while uncovered_pairs:
            best_test = None
            best_score = -1.0
            best_new_pairs = set()

            for i in range(self.numTests):
                if i in chosen_tests or i in excluded_tests:
                    continue

                new_pairs = self.test_to_pairs[i] & uncovered_pairs
                if not new_pairs:
                    continue

                score = score_fn(i, new_pairs)

                if score > best_score:
                    best_test = i
                    best_score = score
                    best_new_pairs = new_pairs

            if best_test is None:
                return None, None

            chosen_tests.add(best_test)
            uncovered_pairs -= best_new_pairs

        chosen_tests = self.remove_redundant_tests(chosen_tests)

        solution_vector = [1 if i in chosen_tests else 0 for i in range(self.numTests)]
        objective_value = sum(self.costOfTest[i] for i in chosen_tests)

        return solution_vector, float(objective_value)

    def greedy_initial_solution(self):
        def greedy_score(i, new_pairs):
            return len(new_pairs) / float(self.costOfTest[i])

        return self.build_feasible_solution(
            initial_chosen_tests=set(),
            excluded_tests=set(),
            score_fn=greedy_score,
        )
    
    def heuristic_solution_from_lp(self, x_values, fixed_zero, fixed_one):
        def lp_guided_score(i, new_pairs):
            return x_values[i] + len(new_pairs) / float(self.costOfTest[i])

        return self.build_feasible_solution(
            initial_chosen_tests=set(fixed_one),
            excluded_tests=set(fixed_zero),
            score_fn=lp_guided_score,
        )

    def forced_one_cost(self, fixed_one):
        return sum(self.costOfTest[i] for i in fixed_one)

    def propagate_forced_choices(self, fixed_zero, fixed_one):
        fixed_zero = set(fixed_zero)
        fixed_one = set(fixed_one)

        changed = True
        while changed:
            changed = False

            for _, _, distinguishing_tests in self.pair_constraints:
                if any(i in fixed_one for i in distinguishing_tests):
                    continue

                available_tests = [i for i in distinguishing_tests if i not in fixed_zero]

                # exactly one possible remaining test => must take it
                if len(available_tests) == 1:
                    forced_test = available_tests[0]
                    fixed_one.add(forced_test)
                    changed = True

        return fixed_zero, fixed_one


    def branch_and_bound(self, fixed_zero, fixed_one):
        counter = count()
        heap = []

        def push_node(node_fixed_zero, node_fixed_one):
            node_fixed_zero, node_fixed_one = self.propagate_forced_choices(
                node_fixed_zero, node_fixed_one
            )

            if self.forced_one_cost(node_fixed_one) >= self.best_objective:
                return

            result = self.solve_lp_relaxation(node_fixed_zero, node_fixed_one)
            if not result["feasible"]:
                return

            if result["objective_value"] >= self.best_objective:
                return

            heapq.heappush(
                heap,
                (
                    result["objective_value"],
                    next(counter),
                    node_fixed_zero,
                    node_fixed_one,
                    result,
                ),
            )

        # push root node
        push_node(set(fixed_zero), set(fixed_one))

        while heap:
            _, _, current_fixed_zero, current_fixed_one, result = heapq.heappop(heap)

            lp_value = result["objective_value"]
            x_values = result["x_values"]

            if lp_value >= self.best_objective:
                continue

            if result["is_integral"]:
                self.best_objective = lp_value
                self.best_solution = self.extract_integer_solution(x_values)
                continue

            heuristic_solution, heuristic_cost = self.heuristic_solution_from_lp(
                x_values, current_fixed_zero, current_fixed_one
            )
            if heuristic_solution is not None and heuristic_cost < self.best_objective:
                self.best_solution = heuristic_solution
                self.best_objective = heuristic_cost

            branch_var = self.choose_branch_variable(
                x_values, current_fixed_zero, current_fixed_one
            )
            if branch_var is None:
                continue

            # left child
            push_node(current_fixed_zero | {branch_var}, set(current_fixed_one))

            # right child
            push_node(set(current_fixed_zero), current_fixed_one | {branch_var})


    def solve(self):
        """
        Healthcare Analytics Model
        """
        if self.is_infeasible:
            return self.solution, self.objective_value

        greedy_solution, greedy_cost = self.greedy_initial_solution()
        if greedy_solution is None:
            return self.solution, self.objective_value

        self.best_solution = greedy_solution
        self.best_objective = greedy_cost
            
        self.branch_and_bound(set(), set())

        self.solution = self.best_solution
        self.objective_value = self.best_objective

        return self.solution, self.objective_value

    def load_from_file(self, filename: str):
        try:
            with open(filename, "r") as fl:
                self.numTests = int(fl.readline().strip())  # n
                self.numDiseases = int(fl.readline().strip())  # m

                self.costOfTest = np.array([float(i) for i in fl.readline().strip().split()])

                self.A = np.zeros((self.numTests, self.numDiseases))
                for i in range(0, self.numTests):
                    self.A[i, :] = np.array([int(i) for i in fl.readline().strip().split()])

        except Exception as e:
            print(f"Error reading instance file. File format may be incorrect.{e}")
            exit(1)

    def __str__(self):
        out = ""
        out = f"Number of test: {self.numTests}\n"
        out += f"Number of diseases: {self.numDiseases}\n"
        cst_str = " ".join([str(i) for i in self.costOfTest])
        out += f"Cost of tests: {cst_str}\n"
        A_str = "\n".join([" ".join([str(j) for j in self.A[i]]) for i in range(0, self.A.shape[0])])
        out += f"A:\n{A_str}"
        return out
