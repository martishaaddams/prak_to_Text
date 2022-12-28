#pragma once
#include "drawings.h"
std::map<int, int> Cuthill_Mckee(int indexofvert, Eigen::SparseMatrix<int> mat, std::vector<vertex>& coord, std::vector<element>& el, std::vector<edge>& ed);
void quickSort(std::vector<int>& arr, int low, int high, std::vector<vertex> coord);
int partition(std::vector<int>& arr, int low, int high, std::vector<vertex> coord);
bool compare_vertex_deg(int v1, int v2, std::vector<vertex> coord);
void swap(int* a, int* b);
void renumerate(std::vector<vertex>& coord, std::map<int, int> renum);
void renumerate(std::vector<element>& el, std::map<int, int> renum);
void print(std::vector<int> a);
int findIndex(std::vector<std::pair<int, double> > a, int x);