#include <bits/stdc++.h>
#include <vector>
using namespace std;

/*
 * Complete the 'solveByGauss' function below.
 *
 * The function is expected to return a DOUBLE_ARRAY.
 * The function accepts following parameters:
 *  1. INTEGER n
 *  2. 2D_DOUBLE_ARRAY matrix
 */
bool isSolutionExists;
string errorMessage;

vector<double> solveByGauss(int n, vector<vector<double>> matrix) {
    vector<double> solution(n);  // Вектор для хранения решения
    vector<double> residuals(n);  // Вектор для хранения невязок
    
    // Прямой ход метода Гаусса с выбором главного элемента
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        // Поиск строки с максимальным по модулю элементом в текущем столбце
        for (int j = i + 1; j < n; ++j) {
            if (abs(matrix[j][i]) > abs(matrix[maxRow][i])) {
                maxRow = j;
            }
        }
        
        // Проверка на вырожденность матрицы
        if (abs(matrix[maxRow][i]) < 1e-10) {
            isSolutionExists = false;
            errorMessage = "The system has no roots of equations or has an infinite set of them.";
            return solution;
        }
        
        // Перестановка строк
        swap(matrix[i], matrix[maxRow]);
        
        // Прямой ход метода Гаусса
        for (int j = i + 1; j < n; ++j) {
            double factor = matrix[j][i] / matrix[i][i];
            for (int k = i; k <= n; ++k) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
    }
    
    // Обратный ход метода Гаусса
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += matrix[i][j] * solution[j];
        }
        solution[i] = (matrix[i][n] - sum) / matrix[i][i];
    }
    
    // Вычисление невязок
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += matrix[i][j] * solution[j];
        }
        residuals[i] = matrix[i][n] - sum;
    }
    
    // Добавление невязок в конец вектора решения
    solution.insert(solution.end(), residuals.begin(), residuals.end());
    
    return solution;
}
