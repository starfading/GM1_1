//
// Created by starfading on 2020/9/8.
//

#pragma once
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class GM1_1{
public:
    vector<double> GetPredictValue(vector<double> original_arr);

private:
    vector<vector<double>> InverseMatrixBxBTxBTxYN(const vector<double>& original_arr);
    vector<vector<double>> MatrixYn(vector<double> original_arr);
    vector<double> CumulativeArray(vector<double> original_arr);
    vector<double> AverageArray(vector<double> cumulative_arr);
    vector<vector<double>> MatrixB(vector<double> average_arr);
    vector<vector<double>> MatrixBT(vector<vector<double>> matrix_b);
    vector<vector<double>> MatrixBxBT(vector<vector<double>> matrix_b, vector<vector<double>> matrix_bt);
    vector<vector<double>> InverseMatrixBxBT(vector<vector<double>> matrix_b_xbt);
    vector<vector<double>> InverseMatrixBxBTxBT(vector<vector<double>> inverse_matrix_b_xbt,
                                                 vector<vector<double>> matrix_bt);
    vector<vector<double>> InverseMatrixBxBTxBTxYN(vector<vector<double>> inverse_matrix_b_xbt_xbt,
                                                    vector<vector<double>> matrix_yn);

};
