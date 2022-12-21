#include"drawings.h"
#define SIZE 3
double poissonRatio, youngModulus;
//treug setka


using namespace std;
/*class csr
{
	vector<int> data;
	vector<int> ja;//Column indices
	vector<int> ia;//stores cumultive number of non zero elements upt ith row

};*/


class element
{
public:
	int num[SIZE];
	void calculate_stifness_matrix(const Eigen::Matrix3f& D, vector<Eigen::Triplet<float> >& triplets)
	{

	}


};
//draw(element)
class vertex
{
public:
	double x;
	double y;
	double z;
	int degree=0;//amount of connections
	vector<int>connections;
	void print()
	{
		cout << ' ' << x << ' ' << y << ' ' << z << ' ';
	}
	
	//vector<int> in_element_id;
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
map<int,int> Cuthill_Mckee(Eigen::SparseMatrix<int> &mat,int indexofvert,vector<vertex>coord)
{
	Eigen::SparseMatrix<int> res;
	map<int, int> renum;
	renum[indexofvert] = 1;
	


	return renum;
}


int main()
{
	string txt;
	//fptr = fopen("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug.k");
	
	ifstream newfile("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug.k");
	ofstream nf("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug_arr.txt");
	//newfile.open("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug.k")
	vector<vertex> arr_Coord;
	vector<element> arr_elem;
	
	std::vector< Eigen::Triplet<int> > triplets;
	//матрица инцендентности
	while (getline(newfile, txt) )
	{
		if (txt == "*NODE")
		{
			break;
		}
	}
	getline(newfile, txt);
	
	while (getline(newfile, txt) )
	{
		if (txt != "$")
		{


			vertex l;
			string delimiter = " ";
			string token;
			int pos;// = txt.find(delimiter);
			//txt.erase(0, 6);
			while (token== "") {
				pos = txt.find(delimiter);
				token = txt.substr(0, pos);
				txt.erase(0, pos + delimiter.length());
			}
			
			int num = stoi(token);
			if (arr_Coord.size() <= num)
			{
				arr_Coord.resize(num + 1);
			}
			token = "";
			while (token == "") {
				pos = txt.find(delimiter);
				token = txt.substr(0, pos);
				txt.erase(0, pos + delimiter.length());
			}
			l.x = stod(token);
			token = "";
			while (token == "") {
				pos = txt.find(delimiter);
				token = txt.substr(0, pos);
				txt.erase(0, pos + delimiter.length());
			}
			l.y = stod(token);
			token = "";
			while (token == "") {
				pos = txt.find(delimiter);
				token = txt.substr(0, pos);
				txt.erase(0, pos + delimiter.length());
			}
			l.z = stod(token);
			arr_Coord[num] = l;
		}
		else break;


	}
	for (int i = 1; i < arr_Coord.size(); i++) {
		nf << arr_Coord[i].x << " " << arr_Coord[i].y << " " << arr_Coord[i].z << endl;

	}nf.close();
	nf.open("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug_arr_nodes.txt");
	getline(newfile, txt);
	getline(newfile, txt);
	getline(newfile, txt);
	
	
	while (getline(newfile, txt))
	{
		if (txt != "*END")
		{
			int n;
			string delimiter = " ";
			string token;
			short int pos;
			while (token == "") {
				pos = txt.find(delimiter);
				token = txt.substr(0, pos);
				txt.erase(0, pos + delimiter.length());
			}
			token = "";
			while (token == "") {
				pos = txt.find(delimiter);
				token = txt.substr(0, pos);
				txt.erase(0, pos + delimiter.length());
			}
			token = "";
			while (token == "") {
				pos = txt.find(delimiter);
				token = txt.substr(0, pos);
				txt.erase(0, pos + delimiter.length());
			}
			n= stoi(token);
			arr_elem.resize(arr_elem.size() + 1);
			nf << n << ' ';
			arr_elem[arr_elem.size() - 1].num[0] = n;
			//arr_Coord[n].degree++;
			//arr_Coord[n].in_element_id[arr_Coord[n].in_element_id.size() - 1] = arr_elem.size()-1;
			for (int i = 1; i < SIZE ; i++)
			{
				token = "";
				while (token == "") {
					pos = txt.find(delimiter);
					token = txt.substr(0, pos);
					txt.erase(0, pos + delimiter.length());
				}
				n = stoi(token);
				nf << n << ' ';
				arr_elem[arr_elem.size() - 1].num[i] = n;
				//arr_Coord[n].in_element_id.resize(arr_Coord[n].in_element_id.size() + 1);
				//arr_Coord[n].in_element_id[arr_Coord[n].in_element_id.size() - 1] = arr_elem.size() - 1;
				token = "";
			}
			for (int k = 0; k < SIZE; k++)
			{
				int i= arr_elem[arr_elem.size() - 1].num[k];
				int j;
				if (k != SIZE - 1)
					{
						 j= arr_elem[arr_elem.size() - 1].num[k + 1];
						 

					}
						else
					{
						j = arr_elem[arr_elem.size() - 1].num[0];
							//j=
					}
				vector<int>::iterator it=find(arr_Coord[i].connections.begin(), arr_Coord[i].connections.end(), j);

				if (it == arr_Coord[i].connections.end())
				{
					
					arr_Coord[i].connections.resize(arr_Coord[i].connections.size() + 1);
					arr_Coord[i].connections[arr_Coord[i].connections.size() - 1] = j;
					arr_Coord[j].connections.resize(arr_Coord[j].connections.size() + 1);
					arr_Coord[j].connections[arr_Coord[j].connections.size() - 1] = i;

				}
			}
			//bool c=true;
			//for (int k = 0; k < SIZE; k++)
			//{
			//	int i = arr_elem[arr_elem.size() - 1].num[k];
			//	int j;
			//	if (k != SIZE-1)
			//	{
			//		 j= arr_elem[arr_elem.size() - 1].num[k + 1];
			//		 

			//	}
			//	else
			//	{
			//		j = arr_elem[arr_elem.size() - 1].num[0];
			//		//j=
			//	}
			//	
			//	triplets.push_back(Eigen::Triplet<int>(i, j, 1));
			//	triplets.push_back(Eigen::Triplet<int>(j, i, 1));
			//	
			//	
		
			//	}
			
			
			/*while (token == "") {
				pos = txt.find(delimiter);
				token = txt.substr(0, pos);
				txt.erase(0, pos + delimiter.length());
			}
			n = stoi(token);
			nf << n << endl;*/


		}
	}

	nf.close();
	newfile.close();
	

	/*for (int i = 0; i < arr_Coord.size(); i++)
	{
		cout << i;
		arr_Coord[i].inelempr();
	}*/
	// the end of reading the file

	//x=0-the line
	Eigen::SparseMatrix<int> adjacency(arr_Coord.size(), arr_Coord.size());
	int mindegree = arr_Coord.size();
	int idmind;
	for (auto i=0; i< arr_Coord.size(); i++)
	{
		arr_Coord[i].degree = arr_Coord[i].connections.size();
		if (mindegree > arr_Coord[i].degree)
		{
			mindegree = arr_Coord[i].degree;
			idmind = i;
		}
		
		for (int k = 0; k < arr_Coord[i].connections.size(); k++)
		{
			triplets.push_back(Eigen::Triplet<int>(i, arr_Coord[i].connections[k], 1));
		}
	}
	adjacency.setFromTriplets(triplets.begin(), triplets.end());
	//adjacency.insert(0, 0);
	adjacency.makeCompressed();
	cout << "Start" << endl;
	draw_matrix(adjacency,"before");

	cout << "Nodes: " << arr_Coord.size() << endl << "Elememts " << arr_elem.size();
	cout << std::endl;//<< adjacency;

	

	//iterate through adjacency
	/*for (int k = 0; k < adjacency.outerSize(); ++k)
		for (Eigen::SparseMatrix<int>::InnerIterator it(adjacency, k); it; ++it)
		{
			cout<<it.value() << ' ';
			cout<<it.row() << ' ';   // row index
			cout<<it.col() << ' ';  // col index (here it is equal to k)
			cout << endl;
			//cout<<it.index() << ' '; // inner index, here it is equal to it.row()
		}*/

	//elasticity matrix

	Eigen::Matrix3f D;
	D <<
		1.0f, poissonRatio, 0.0f,
		poissonRatio, 1.0, 0.0f,
		0.0f, 0.0f, (1.0f - poissonRatio) / 2.0f;

	D *= youngModulus / (1.0f - pow(poissonRatio, 2.0f));



	/*Eigen::SparseMatrix<int> A(5, 5);

	std::vector< Eigen::Triplet<int> > triplets;

	triplets.push_back(Eigen::Triplet<int>(0, 1, 3));
	triplets.push_back(Eigen::Triplet<int>(1, 0, 22));
	triplets.push_back(Eigen::Triplet<int>(2, 1, 5));
	triplets.push_back(Eigen::Triplet<int>(2, 3, 1));
	triplets.push_back(Eigen::Triplet<int>(4, 2, 14));
	triplets.push_back(Eigen::Triplet<int>(4, 4, 8));

	A.setFromTriplets(triplets.begin(), triplets.end());

	A.insert(0, 0);
	std::cout << A;

	A.makeCompressed();

	std::cout << std::endl << A;*/
	


	return 0;

	


}