#include <bits/stdc++.h>
using namespace std;

auto linearApproximation(const vector<double> &x, const vector<double> &y)
    -> double {
    auto n = x.size();

    if (n < 2) {
        return 0.0; // Недостаточно точек для аппроксимации
    }

    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    double denom = n * sumX2 - sumX * sumX;
    if (abs(denom) < 1e-10) {
        return 0.0; // Деление на ноль, точки лежат на вертикальной прямой
    }

    double a = (n * sumXY - sumX * sumY) / denom;
    double b = (sumY * sumX2 - sumX * sumXY) / denom;

    double maxSqDev = 0.0;
    for (int i = 0; i < n; i++) {
        double approxY = a * x[i] + b;
        double sqDev = pow(y[i] - approxY, 2);
        maxSqDev = max(maxSqDev, sqDev);
    }

    return maxSqDev;
}

auto main() -> int {
    int n;
    cin >> n;

    vector<double> x(n), y(n);
    for (auto &xi : x) {
        cin >> xi;
    }
    for (auto &yi : y) {
        cin >> yi;
    }

    double maxSqDeviation = linearApproximation(x, y);
    cout << fixed << setprecision(10) << maxSqDeviation << endl;
}
