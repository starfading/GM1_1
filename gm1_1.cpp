//
// Created by starfading on 2020/9/8.
//

#include "gm1_1.h"

vector<double> GM1_1::GetPredictValue(vector<double> original_arr) {
    int n = original_arr.size();
    if(n <= 1) {
        cerr << "original array is invalid" << '\n';
        exit(0);
    }

    // 获得预测参数：a,u
    vector<vector<double>> inverse_matrix_b_xbt_xbt_xyn = InverseMatrixBxBTxBTxYN(original_arr);
    double a = inverse_matrix_b_xbt_xbt_xyn[0][0];
    double u = inverse_matrix_b_xbt_xbt_xyn[1][0];
    double au = 0;
    if(a != 0) au = u / a;

    // 获得预测序列
    int m = 5;
    vector<double> predict_result(n + m);
    predict_result[0] = original_arr[0];
    for(int i = 1; i < n + m; ++i){
        predict_result[i] = (original_arr[0] - au) * exp(-a * i) * (1 - exp(a));
    }
    cout << "predict  array: ";
    for(auto &v: predict_result) cout << v << '\t';

    // 预测误差评定
    double avg = 0;
    for(auto &v: original_arr) avg += v;
    avg /= n;
    double er_max = 0;
    for(int i = 0; i < n; ++i){
        double tmp = abs(original_arr[i] - predict_result[i]) / avg;
        if(tmp > er_max) er_max = tmp;
    }
    if(er_max <= 0.15) cout << '\n' << "great predict, " << "error ratio is " << er_max << '\n';
    else cout << '\n' << "bad predict, " << "error ratio is " << er_max << '\n';

    return predict_result;
}

vector<vector<double>> GM1_1::InverseMatrixBxBTxBTxYN(const vector<double>& original_arr) {
    vector<vector<double>> matrix_yn = MatrixYn(original_arr);
    vector<double> cumulative_arr = CumulativeArray(original_arr);
    vector<double> average_arr = AverageArray(cumulative_arr);
    vector<vector<double>> matrix_b = MatrixB(average_arr);
    vector<vector<double>> matrix_bt = MatrixBT(matrix_b);
    vector<vector<double>> matrix_b_xbt = MatrixBxBT(matrix_b, matrix_bt);
    vector<vector<double>> inverse_matrix_b_xbt = InverseMatrixBxBT(matrix_b_xbt);
    vector<vector<double>> inverse_matrix_b_xbt_xbt = InverseMatrixBxBTxBT(inverse_matrix_b_xbt, matrix_bt);
    vector<vector<double>> inverse_matrix_b_xbt_xbt_yn = InverseMatrixBxBTxBTxYN(inverse_matrix_b_xbt_xbt, matrix_yn);
    return inverse_matrix_b_xbt_xbt_yn;
}

vector<vector<double>> GM1_1::MatrixYn(vector<double> original_arr) {
    vector<vector<double>> matrix_yn(vector<vector<double>>(original_arr.size()-1, vector<double>(1)));
    for (int i = 0; i < matrix_yn.size(); ++i) {
        matrix_yn[i][0] = original_arr[i + 1];
    }
    return matrix_yn;
}

vector<double> GM1_1::CumulativeArray(vector<double> original_arr) {
    vector<double> cumulative_arr(original_arr.size());
    double sum = 0;
    for (int i = 0; i < cumulative_arr.size(); ++i) {
        sum += original_arr[i];
        cumulative_arr[i] = sum;
    }
    return cumulative_arr;
}

vector<double> GM1_1::AverageArray(vector<double> cumulative_arr) {
    vector<double> average_arr(cumulative_arr.size()-1);
    for (int i = 0; i < average_arr.size(); ++i) {
        average_arr[i] = (cumulative_arr[i] + cumulative_arr[i + 1]) / 2.0;
    }
    return average_arr;
}

vector<vector<double>> GM1_1::MatrixB(vector<double> average_arr) {
    vector<vector<double>> matrix_b(vector<vector<double>>(average_arr.size(), vector<double>(2)));
    for (int i = 0; i < matrix_b.size(); ++i) {
        matrix_b[i][0] = (-average_arr[i]);
        matrix_b[i][1] = 1.0;
    }
    return matrix_b;
}

vector<vector<double>> GM1_1::MatrixBT(vector<vector<double>> matrix_b) {
    vector<vector<double>> matrix_bt(vector<vector<double>>(2, vector<double>(matrix_b.size())));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < matrix_b.size(); ++j) {
            matrix_bt[i][j] = matrix_b[j][i];
        }
    }
    return matrix_bt;
}

vector<vector<double>> GM1_1::MatrixBxBT(vector<vector<double>> matrix_b, vector<vector<double>> matrix_bt) {
    vector<vector<double>> matrix_b_xbt(vector<vector<double>>(2, vector<double>(2)));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < matrix_b.size(); ++k) {
                matrix_b_xbt[i][j] += matrix_bt[i][k] * matrix_b[k][j];
            }
        }
    }
    return matrix_b_xbt;
}

vector<vector<double>> GM1_1::InverseMatrixBxBT(vector<vector<double>> matrix_b_xbt) {
    if (matrix_b_xbt.empty()) return matrix_b_xbt;
    vector<vector<double>> inverse_matrix_b_xbt(vector<vector<double>>(2, vector<double>(2)));
    double tmp = matrix_b_xbt[0][0] * matrix_b_xbt[1][1] - matrix_b_xbt[0][1] * matrix_b_xbt[1][0];
    inverse_matrix_b_xbt[0][0] = 1.0 / tmp * matrix_b_xbt[1][1];
    inverse_matrix_b_xbt[0][1] = 1.0 / tmp * -matrix_b_xbt[0][1];
    inverse_matrix_b_xbt[1][0] = 1.0 / tmp * -matrix_b_xbt[1][0];
    inverse_matrix_b_xbt[1][1] = 1.0 / tmp * matrix_b_xbt[0][0];
    return inverse_matrix_b_xbt;
}

vector<vector<double>> GM1_1::InverseMatrixBxBTxBT(vector<vector<double>> inverse_matrix_b_xbt, vector<vector<double>> matrix_bt) {
    if (inverse_matrix_b_xbt.empty() || matrix_bt.empty()) return inverse_matrix_b_xbt;
    vector<vector<double>> inverse_matrix_b_xbt_xbt(vector<vector<double>>(2, vector<double>(matrix_bt[0].size())));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < matrix_bt[0].size(); ++j) {
            for (int k = 0; k < 2; ++k) {
                inverse_matrix_b_xbt_xbt[i][j] += inverse_matrix_b_xbt[i][k] * matrix_bt[k][j];
            }
        }
    }
    return inverse_matrix_b_xbt_xbt;
}

vector<vector<double>> GM1_1::InverseMatrixBxBTxBTxYN(vector<vector<double>> inverse_matrix_b_xbt_xbt, vector<vector<double>> matrix_yn) {
    if (inverse_matrix_b_xbt_xbt.empty() || matrix_yn.empty()) return inverse_matrix_b_xbt_xbt;
    vector<vector<double>> inverse_matrix_b_xbt_xbt_xyn(vector<vector<double>>(2, vector<double>(1)));
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 1; ++j) {
            for (int k = 0; k < matrix_yn.size(); ++k) {
                inverse_matrix_b_xbt_xbt_xyn[i][j] += inverse_matrix_b_xbt_xbt[i][k] * matrix_yn[k][j];
            }
        }
    }
    return inverse_matrix_b_xbt_xbt_xyn;
}
