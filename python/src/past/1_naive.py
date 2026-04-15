from typing import Any


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
        
        solver = pywraplp.Solver.CreateSolver("GLOP")
        
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
        pair_constraints = self.get_distinguishing_test_sets()
        for j, k, distinguishing_tests in pair_constraints:

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
        for i in range(self.numTests):
            if i in fixed_zero or i in fixed_one:
                continue
            if abs(x_values[i] - round(x_values[i])) > eps:
                return i
        return None

    def extract_integer_solution(self, x_values):
        return [round(v) for v in x_values]

    def branch_and_bound(self, fixed_zero, fixed_one):
        result = self.solve_lp_relaxation(fixed_zero, fixed_one)
        
        if not result["feasible"]:
            return
        
        lp_value = result["objective_value"]
        x_values = result["x_values"]

        if lp_value >= self.best_objective:
            return
        
        if result["is_integral"]:
            self.best_objective = lp_value
            self.best_solution = self.extract_integer_solution(x_values)
            return
        
       
        i = self.choose_branch_variable(x_values, fixed_zero, fixed_one)
        # should never be None
        if i is None:
            print("ERROR: should never reach here")
            return

        self.branch_and_bound(fixed_zero | {i}, fixed_one)
        self.branch_and_bound(fixed_zero, fixed_one | {i})
        
    def solve(self):
        """
        Healthcare Analytics Model
        """

        result = self.branch_and_bound(set(), set())

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
