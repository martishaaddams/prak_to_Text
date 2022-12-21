#pragma once
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<map>
#include<cmath>
#include<Eigen/Dense>
#include<Eigen/Sparse>
int draw_matrix(Eigen::SparseMatrix<int> mat,std::string fn="none");
