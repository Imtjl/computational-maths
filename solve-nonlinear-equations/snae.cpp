#include <vector>
using namespace std;

const double STEP_SIZE = 1e-8;
const double TOLERANCE = 1e-10;
const double ROUND_FACTOR = 1e5;
const int MAX_ITERATIONS = 10000;

vector<double> solveByGauss(int n, vector<vector<double>> &matrix) {
    vector<double> solution(n, 0.0);

    // Gauss elimination with pivoting
    for (int i = 0; i < n; ++i) {
        // Finding the row with the maximum element in the current column.
        int maxRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(matrix[j][i]) > abs(matrix[maxRow][i])) {
                maxRow = j;
            }
        }

        if (abs(matrix[maxRow][i]) < 1e-10) { // Check for zero pivot element.
            return vector<double>(); // Return an empty vector to indicate no
                                     // solution.
        }

        // Swapping the maximum row with the current row.
        swap(matrix[i], matrix[maxRow]);

        // Forward elimination
        for (int j = i + 1; j < n; ++j) {
            double factor = matrix[j][i] / matrix[i][i];
            for (int k = i; k <= n; ++k) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
    }

    // Back Substitution
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += matrix[i][j] * solution[j];
        }
        solution[i] = (matrix[i][n] - sum) / matrix[i][i];
    }
    return solution;
}

vector<vector<double>> computeJacobian(const funcvec_t &functions,
                                       const vector<double> &values) {
    int n = values.size();
    vector<vector<double>> jacobian(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        vector<double> derivatives(n);
        double function_value = (*functions[i])(values);
        for (int j = 0; j < n; j++) {
            vector<double> perturbed_values = values;
            perturbed_values[j] += STEP_SIZE;
            double perturbed_function_value = (*functions[i])(perturbed_values);
            derivatives[j] =
                (perturbed_function_value - function_value) / STEP_SIZE;
        }
        jacobian[i] = derivatives;
    }
    return jacobian;
}

vector<double> solve_by_newton(int system_id, int num_unknowns,
                               vector<double> initial_values) {
    funcvec_t functions = get_functions(system_id);
    if (num_unknowns != initial_values.size() ||
        functions.size() != initial_values.size()) {
        return vector<double>();
    }
    vector<double> current_values = initial_values;
    double max_difference;
    int iterations = 0;
    double damping_factor = 1.0;
    do {
        vector<double> function_values(functions.size());
        for (int i = 0; i < functions.size(); i++) {
            function_values[i] = (*functions[i])(current_values);
        }
        double max_function_value =
            *max_element(function_values.begin(), function_values.end(),
                         [](double a, double b) { return abs(a) < abs(b); });
        if (max_function_value < TOLERANCE) {
            // Function values are sufficiently small, consider it as a solution
            break;
        }

        vector<vector<double>> jacobian_matrix =
            computeJacobian(functions, current_values);
        for (int i = 0; i < functions.size(); i++) {
            function_values[i] = -function_values[i];
        }
        for (int i = 0; i < jacobian_matrix.size(); ++i) {
            jacobian_matrix[i].push_back(function_values[i]);
        }
        vector<double> delta_values =
            solveByGauss(jacobian_matrix.size(), jacobian_matrix);
        if (delta_values.empty()) {
            // Singular Jacobian matrix, perturb the current approximation
            // slightly
            for (int i = 0; i < current_values.size(); i++) {
                current_values[i] += STEP_SIZE;
            }
            continue;
        }
        for (int i = 0; i < current_values.size(); i++) {
            current_values[i] += damping_factor * delta_values[i];
        }
        max_difference =
            *max_element(delta_values.begin(), delta_values.end(),
                         [](double a, double b) { return abs(a) < abs(b); });
        iterations++;

        // Adjust the damping factor based on the progress
        if (max_difference < TOLERANCE) {
            damping_factor *= 2.0;
        } else {
            damping_factor *= 0.5;
        }
    } while (max_difference >= TOLERANCE && iterations < MAX_ITERATIONS);

    if (iterations == MAX_ITERATIONS) {
        // Maximum iterations reached without convergence
        return vector<double>();
    }

    for (int i = 0; i < current_values.size(); i++) {
        current_values[i] =
            round(current_values[i] * ROUND_FACTOR) / ROUND_FACTOR;
    }
    return current_values;
}
