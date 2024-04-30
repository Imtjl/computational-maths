#include <math.h>
#include <string>
using namespace std;

string error_message = "";
bool has_discontinuity = false;

double first_function(double x) { return 1 / x; }

double second_function(double x) { return sin(x) / x; }

double third_function(double x) { return x * x + 2; }

double fourth_function(double x) { return 2 * x + 2; }

double five_function(double x) { return log(x); }

/*
 * How to use this function:
 *    fn_t& f = get_function(4);
 *    f(0.0001);
 */
fn_t &get_function(int n) {
    switch (n) {
    case 1:
        return first_function;
    case 2:
        return second_function;
    case 3:
        return third_function;
    case 4:
        return fourth_function;
    case 5:
        return five_function;
    default:
        throw std::invalid_argument("Invalid function number");
    }
}

/*
 * Complete the 'calculate_integral' function below.
 *
 * The function is expected to return a DOUBLE.
 * The function accepts following parameters:
 *  1. DOUBLE a
 *  2. DOUBLE b
 *  3. INTEGER f
 *  4. DOUBLE epsilon
 */

double calculate_integral(double lower_bound, double upper_bound, int f,
                          double epsilon) {
    fn_t &function = get_function(f);

    // Check for discontinuity or undefined region
    if ((f == 1 && (lower_bound <= 0 || upper_bound <= 0)) ||
        (f == 5 && (lower_bound <= 0 || upper_bound <= 0))) {
        error_message = "Integrated function has discontinuity or does not "
                        "defined in current interval";
        has_discontinuity = true;
        return 0.0;
    }

    int num_trapezoids = 1;
    double integral = 0.0;
    double prev_integral = 0.0;

    do {
        prev_integral = integral;
        integral = 0.0;
        double trapezoid_width = (upper_bound - lower_bound) / num_trapezoids;

        for (int i = 0; i < num_trapezoids; i++) {
            double x1 = lower_bound + i * trapezoid_width;
            double x2 = lower_bound + (i + 1) * trapezoid_width;
            integral += 0.5 * (function(x1) + function(x2)) * trapezoid_width;
        }

        num_trapezoids *= 2;
    } while (std::abs(integral - prev_integral) > epsilon);

    return integral;
}
