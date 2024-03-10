#include <bits/stdc++.h>
#include <vector>
using namespace std;

/*
 * Complete the 'approximate_linear_least_squares' function below.
 *
 * The function is expected to return a DOUBLE.
 * The function accepts following parameters:
 *  1. DOUBLE_ARRAY x_axis
 *  2. DOUBLE_ARRAY y_axis
 */
double approximate_linear_least_squares(vector<double> xVals,
                                        vector<double> yVals) {
    int n = xVals.size();

    // Для аппроксимации должно быть минимум две точки
    if (n < 2)
        return 0.0;

    // Вычисляю сумму для коэффициентов
    double xSum = accumulate(xVals.begin(), xVals.end(), 0.0);
    double ySum = accumulate(yVals.begin(), yVals.end(), 0.0);
    double xySum =
        inner_product(xVals.begin(), xVals.end(), yVals.begin(), 0.0);
    double xSqSum =
        inner_product(xVals.begin(), xVals.end(), xVals.begin(), 0.0);

    // Вычисляю знаменатель для коэффициентов
    double denominator = n * xSqSum - xSum * xSum;

    // Если знаменатель близок к нулю, то точки лежат на вертикальной прямой
    if (abs(denominator) < 1e-10)
        return 0.0;

    // Вычисление коэффициентов для аппроксимирующей прямой
    double aCoeff = (n * xySum - xSum * ySum) / denominator;
    double bCoeff = (ySum * xSqSum - xSum * xySum) / denominator;

    // Нахожу максимальное квадратичное отклонение
    double maxSqDev = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double approximateY = aCoeff * xVals[i] + bCoeff;
        double sqDev = pow(yVals[i] - approximateY, 2);
        maxSqDev = max(maxSqDev, sqDev);
    }

    return maxSqDev;
}
