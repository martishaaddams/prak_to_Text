#pragma once
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<queue>
#include<map>
#include<cmath>
#include<complex>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#define SIZE 3
int draw_matrix(Eigen::SparseMatrix<int> mat,std::string fn="none");

class vertex
{
public:
	double x;
	double y;
	double z;
	std::vector<int> in_element_id;
	int degree = 0;//amount of connections
	std::vector<int>connections;
	void print()
	{
		std::cout << ' ' << x << ' ' << y << ' ' << z << ' ';
	}
	void renum(std::map<int, int> renum)
	{
		for (auto i = 0; i < connections.size(); i++)
		{
			if (renum.find(connections[i]) != renum.end())
				connections[i] = renum[connections[i]];
			else
				std::cout << "Something is wrong in remum in vertex" << std::endl;
		}
	}
	void renum(std::vector<int> renum)
	{
		for (auto i = 0; i < connections.size(); i++)
		{
			connections[i] = renum[connections[i]];
		}
	}


	
	/*void inelempr()
	{
		cout << endl;
		for (int j = 0; j < in_element_id.size(); j++)
		{

			cout << in_element_id[j] << ' ';
		}
		cout << endl;
	}*/
};
class element
{
public:
	int num[SIZE];
	/*void calculate_stifness_matrix(const Eigen::Matrix3f& D, std::vector<Eigen::Triplet<float> >& triplets)
	{

	}*/
	void renum(std::map<int, int> renum)
	{
		for (auto i = 0; i < SIZE; i++)
		{
			if (renum.find(num[i]) != renum.end())
				num[i] = renum[num[i]];
			else
				std::cout << "Something is wrong in remum in element" << std::endl;
			
		}
	}
	void renum(std::vector<int> renum)
	{
		for (auto i = 0; i < SIZE; i++)
		{
			num[i] = renum[num[i]];
		}
	}
	//onlyfortriangle
	std::vector<Eigen::Triplet<float> > CalculateStiffnessMatrix(const Eigen::Matrix3f& D, std::vector<Eigen::Triplet<float> > triplets, std::vector<vertex> coord)
	{
		Eigen::Vector3f x, y;
		/*x << coord[num[0]].x, coord[num[1]].x, coord[num[2]].x;
		y << coord[num[0]].y, coord[num[1]].y, coord[num[2]].y;*/
		for (int i = 0; i < SIZE; i++) {
			x[i] = coord[num[i]].x;
			y[i] = coord[num[i]].y;
		}
		Eigen::Matrix3f C;
		C << Eigen::Vector3f(1.0f, 1.0f, 1.0f), x, y;
		//std::cout << C.determinant();
		Eigen::Matrix3f IC = C.inverse();
		Eigen::Matrix<float, 3, 6> B;
		for (int i = 0; i < 3; i++)
		{
			B(0, 2 * i + 0) = IC(1, i);
			B(0, 2 * i + 1) = 0.0f;
			B(1, 2 * i + 0) = 0.0f;
			B(1, 2 * i + 1) = IC(2, i);
			B(2, 2 * i + 0) = IC(2, i);
			B(2, 2 * i + 1) = IC(1, i);
		}
		//std::cout << B;
		Eigen::Matrix<float, 6, 6> K = B.transpose() * D * B * C.determinant() / 2.0f;
		//std::cout <<  D * B  << std::endl;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				Eigen::Triplet<float> trplt11(2 * num[i] + 0, 2 * num[j] + 0, K(2 * i + 0, 2 * j + 0));
				Eigen::Triplet<float> trplt12(2 * num[i] + 0, 2 * num[j] + 1, K(2 * i + 0, 2 * j + 1));
				Eigen::Triplet<float> trplt21(2 * num[i] + 1, 2 * num[j] + 0, K(2 * i + 1, 2 * j + 0));
				Eigen::Triplet<float> trplt22(2 * num[i] + 1, 2 * num[j] + 1, K(2 * i + 1, 2 * j + 1));

				triplets.push_back(trplt11);
				triplets.push_back(trplt12);
				triplets.push_back(trplt21);
				triplets.push_back(trplt22);
			}
		}
		/*for (int i = 0; i < triplets.size(); i++)
		{
			std::cout << triplets[i].row() << " " << triplets[i].col() << " " << triplets[i].value() << std::endl;
		}*/
		return triplets;
	}
	
	


};

class edge
{
public:
	int start;
	int end;
	int in_element;
	void renumerate(std::vector<int> renum)
	{
		start = renum[start];
		end = renum[end];
	}
	
};
int draw_edges(std::vector<edge> ed, std::vector<vertex> coord, std::string fn="none");
int draw_coord(std::vector<vertex> coord);
int draw_graph(std::vector<float> x, std::vector<float> y);