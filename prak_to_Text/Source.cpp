//#include"drawings.h"
#include"somealgo.h"
double poissonRatio=0.25;
double youngModulus = 2 * pow(10, 11);
double PP = pow(10, 8);
//treug setka
//E = 2 * 10 * *11
//nu = 0.25
//P = 10 * *8

using namespace std;
/*class csr
{
	vector<int> data;
	vector<int> ja;//Column indices
	vector<int> ia;//stores cumultive number of non zero elements upt ith row

};*/
/*
void print(vector<int> a)
{
	for (auto i = 0; i < a.size(); i++)
	{
		cout << a[i] << ' ';
	}
	cout << endl;
}
void swap(int* a, int *b)
{
	int t = *a;
	*a = *b;
	*b = t;
}
bool compare_vertex_deg(int v1,int v2,vector<vertex> coord)
{
	return(coord[v1].degree < coord[v2].degree);
	//return(v1.x < v2.x);
}
/*void merge(vector<int>& array, int const left, int const mid,
	int const right, vector<vertex> coord)
{
	auto const subArrayOne = mid - left + 1;
	auto const subArrayTwo = right - mid;

	// Create temp arrays
	//auto* leftArray = new int[subArrayOne],
	//	* rightArray = new int[subArrayTwo];
	vector<int> leftArray;
	leftArray.resize(subArrayOne);
	vector<int> rightArray;
	rightArray.resize(subArrayOne);
	// Copy data to temp arrays leftArray[] and rightArray[]
	for (auto i = 0; i < subArrayOne; i++)
		leftArray[i] = array[left + i];
	for (auto j = 0; j < subArrayTwo; j++)
		rightArray[j] = array[mid + 1 + j];

	auto indexOfSubArrayOne
		= 0, // Initial index of first sub-array
		indexOfSubArrayTwo
		= 0; // Initial index of second sub-array
	int indexOfMergedArray
		= left; // Initial index of merged array

	// Merge the temp arrays back into array[left..right]
	while (indexOfSubArrayOne < subArrayOne
		&& indexOfSubArrayTwo < subArrayTwo) {
		if (!compare_vertex_deg(rightArray[indexOfSubArrayTwo], leftArray[indexOfSubArrayOne],coord))
			 {
			array[indexOfMergedArray]
				= leftArray[indexOfSubArrayOne];
			indexOfSubArrayOne++;
		}
		else {
			array[indexOfMergedArray]
				= rightArray[indexOfSubArrayTwo];
			indexOfSubArrayTwo++;
		}
		indexOfMergedArray++;
	}
	// Copy the remaining elements of
	// left[], if there are any
	while (indexOfSubArrayOne < subArrayOne) {
		array[indexOfMergedArray]
			= leftArray[indexOfSubArrayOne];
		indexOfSubArrayOne++;
		indexOfMergedArray++;
	}
	// Copy the remaining elements of
	// right[], if there are any
	while (indexOfSubArrayTwo < subArrayTwo) {
		array[indexOfMergedArray]
			= rightArray[indexOfSubArrayTwo];
		indexOfSubArrayTwo++;
		indexOfMergedArray++;
	}
	
}
void merge_sort_connections(int l, int r, vector<int>& con,vector<vertex> coord)
{
	if (r >= l)
	{
		
		return ;
	}
	
		int m = (r + l) / 2;
		merge_sort_connections(l, m,con,coord);
		merge_sort_connections(m+1,r,con,coord);
		merge(con, l, m, r,coord);

	
	return ;

}
int partition(vector<int>& arr, int low, int high, vector<vertex> coord)
{
	int pivot = coord[arr[high]].degree; // pivot
	int i
		= (low
			- 1); // Index of smaller element and indicates
				  // the right position of pivot found so far

	for (int j = low; j <= high - 1; j++) {
		// If current element is smaller than the pivot
		if (coord[arr[j]].degree < pivot) {
			i++; // increment index of smaller element
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);
}
void quickSort(vector<int>& arr, int low, int high, vector<vertex> coord)
{
	if (low < high) {
		/* pi is partitioning index, arr[p] is now
		at right place 
		int pi = partition(arr, low, high,coord);

		// Separately sort elements before
		// partition and after partition
		quickSort(arr, low, pi - 1,coord);
		quickSort(arr, pi + 1, high,coord);
	}
}*/
Eigen::SparseMatrix<int> build_matrix(vector<vertex>& coord)
{
	Eigen::SparseMatrix<int> res(coord.size(),coord.size());
	vector<Eigen::Triplet<int>> triplets;
	for (auto i = 0; i < coord.size(); i++)
	{
		for (int k = 0; k < coord[i].connections.size(); k++)
		{
			triplets.push_back(Eigen::Triplet<int>(i, coord[i].connections[k], 1));
		}
	}
	res.setFromTriplets(triplets.begin(), triplets.end());
	//adjacency.insert(0, 0);
	res.makeCompressed();
	return res;
}
/*void renumerate(vector<vertex>& coord, map<int, int> renum)
{
	vertex temp;
	for (int i = 0; i < coord.size(); i++)
	{
		temp = coord[i];
		coord[i] = coord[renum[i]];
		coord[i] = temp;
	}
	for (int i = 0; i < coord.size(); i++)
	{
		coord[i].renum(renum);
	}
}
void renumerate(vector<element>& el, map<int, int> renum)
{
	for (int i = 0; i < el.size(); i++)
	{
		el[i].renum(renum);
	}
}*/
//void sort(std::vector int)
/*map<int, int> Cuthill_Mckee(int indexofvert, Eigen::SparseMatrix<int> mat, vector<vertex>& coord, vector<element>& el)
{
	Eigen::SparseMatrix<int> res(mat.rows(), mat.cols());

	map<int, int> renum;
	queue<int> qu;
	vector<int> renumx;
	////renum[indexofvert] = 1;
	qu.push(indexofvert);

	int k = 1;
	while (!qu.empty())
	{
		if (renum.find(qu.front()) == renum.end())
		{
			/*print(coord[qu.front()].connections);
			for (int i = 0; i < coord[qu.front()].connections.size(); i++)
			{
				cout << coord[coord[qu.front()].connections[i]].degree << ' ';
			}
			cout << endl;
			//auto it = coord[qu.front()].connections.begin()
			//merge_sort_connections(0,coord[ qu.front()].connections.size(), coord[qu.front()].connections,coord);
			//swap(&coord[qu.front()].connections[0], &coord[qu.front()].connections[coord[qu.front()].connections.size() - 1]);
			quickSort(coord[qu.front()].connections, 0, coord[qu.front()].connections.size() - 1,coord);
			/*print(coord[qu.front()].connections);
			for (int i = 0; i < coord[qu.front()].connections.size();i++)
			{
				cout << coord[coord[qu.front()].connections[i]].degree << ' ';
			}
			cout << endl;
			print(coord[qu.front()].connections);
			for (int i = 0; i < coord[qu.front()].connections.size(); i++)
			{
				cout << coord[coord[qu.front()].connections[i]].degree << ' ';
			}
			cout << endl;
			for (auto i = 0; i < coord[qu.front()].connections.size(); i++)
			{
				qu.push(coord[qu.front()].connections[i]);
			}

			renum[qu.front()] = k;
			renumx[k] = qu.front();
			k++;
			cout << qu.front() << endl;
		}
		qu.pop();
	}
	
	renumerate(coord, renum);
	renumerate(el, renum);
	return renum;
}*/
int closest(double x, double y, vector<vertex> c)
{
	double r;
	double rmin=-1;
	int index;
	for (int i = 1; i < c.size(); i++)
	{
		
		r = sqrt(pow(x - c[i].x, 2) + pow(y - c[i].y, 2));
		if (rmin == -1 || r < rmin)
		{
			rmin = r;
			index = i;
		}
		
		
	}
	return index;
}
vector<float> getanalyticalsol(vertex loc, float p, float a)
{
	vector<float> res;
	float r = sqrt(pow(loc.x, 2) + pow(loc.y, 2));
	float thete = acos(loc.x / r);
	//res.resize(3);
	/*
	res[0] = p * (1 - pow(a, 2) / pow(r, 2)) / 2 + p * (1 - 4 * pow(a, 2) / pow(r, 2) + 3 * pow(a, 4) / pow(r, 4)) * cos(2 * thete) / 2;
	res[1] = p * (1 + pow(a, 2) / pow(r, 2)) / 2 - p * (1 + 3 * pow(a, 4) / pow(r, 4)) * cos(2 * thete) / 2;
	res[2] = p * (1 + 2 * pow(a, 2) / pow(r, 2) - 3 * pow(a, 4) / pow(r, 4)) * sin(thete) / 2;

	dphi  = np.stack(1/4*p0*(1+2*a**2/(r*np.exp(theta*1j))**2),axis=0);
        d2phi = -p0*a**2/(r*np.exp(1j*theta))**3;
        dpsi  = 1/2*p0*(1+a*(a/(r*np.exp(1j*theta))**2 + 3*a**3/(r*np.exp(1j*theta))**4));

        p_theta = (dphi + np.conj(dphi)) + (np.exp(1j*theta)**2*(np.conj(r*np.exp(1j*theta))*d2phi + dpsi)).real;
        p_r     = (dphi + np.conj(dphi)) - (np.exp(1j*theta)**2*(np.conj(r*np.exp(1j*theta))*d2phi + dpsi)).real;
        p_r_theta = ((np.exp(1j*theta))**2*((np.conj(r*np.exp(1j*theta)))*d2phi + dpsi)).imag;
        return p_r,p_theta,p_r_theta
	*/
	std::complex<float> rphi(cos(thete), sin(thete));// = polar(1.0f, thete);e^iphi
	std::complex<float> dphi = p * (1.0f + 2.0f * complex<float>( pow(a, 2),0) / complex<float>(pow((rphi), 2))) / 4.0f;
	std::complex<float> d2phi = -p * complex<float>(pow(a, 2), 0) / complex<float>(pow((r * rphi), 3));
	std::complex<float> dpsi = 1 / 2 * p * (1.0f + a * (a / ((r * rphi) * (r * rphi)) + 3.0f * complex<float>(pow(a, 3)) / ((r * rphi) * (r * rphi) * (r * rphi) * (r * rphi))));
	
	std::complex<float> p_theta = (dphi + conj(dphi)) + real(rphi * rphi * (conj(r * rphi) * d2phi + dpsi));
	std::complex<float> p_r = (dphi + conj(dphi)) - real(rphi * rphi * (conj(r * rphi)) * d2phi + dpsi);
	std::complex<float> p_r_theta = imag(rphi*rphi*(conj(r*rphi)*d2phi+dpsi));
	res.push_back(real(p_r));
	res.push_back(real(p_theta));
	res.push_back(real(p_r_theta));
	return res;
}
//Eigen::VectorXf interpolation(Eigen)
//vector<float>solution(double xb, double yb, double xe, double ye, double h, Eigen::VectorXf u, Eigen::Matrix3f D)
//{
//	//sigma = dbu
//
//}

/*vector<double> step(double x, double y, vector<vertex>& c, vector<element>& el, Eigen::VectorXf u)//vector<vector<el>> boundaries)
{
	int index = closest(x, y, c);

	
}*/
vector<float> analytical(double xb, double yb,double xe,double ye,double h)
{
	vector<float> res;
	vertex loc;
	double x = xb;
	double y = yb;
	bool c = false;
	if (yb == ye)
	{
		c = true;
	}
	while (x < xe)
	{
		x = x + h;
		if (!c)
		{
			y = (yb - ye) * x / (xe - xb) ;

		}
		loc.x = x;
		loc.y = y;
		res.push_back(getanalyticalsol(loc,PP,0.01)[0]);
			
	}
	return res;
}

vector<edge> getbounds(vector<edge> arredge)
{
	vector<edge> bounds;
	for (int i = 0; i < arredge.size(); i++)
	{
		if (arredge[i].in_element == 1)
		{
			bounds.push_back(arredge[i]);
		}
	}
	return bounds;
}

vector<edge> getleftbounds(vector<edge> b, vector<vertex> c)
{
	vector<edge> leftbounds;
	vector<edge> rightbounds;
	vector<edge> all;
	double minx;
	for (int i = 0; i < b.size(); i++)
	{
		if (i == 0)
		{
			minx = c[b[i].start].x;
		}
		if (c[b[i].start].x == c[b[i].end].x)
		{
			if (c[b[i].start].x < minx)
			{
				minx = c[b[i].start].x;
			}
			all.push_back(b[i]);
		}
	}
	for (int i = 0; i < all.size(); i++)
	{
		if (c[all[i].start].x == minx)
		{
			leftbounds.push_back(all[i]);
		}
	}
	return leftbounds;

}
vector<edge> getrightbounds(vector<edge> b, vector<vertex> c)
{
	vector<edge> leftbounds;
	vector<edge> rightbounds;
	vector<edge> all;
	double maxx;
	for (int i = 0; i < b.size(); i++)
	{
		if (i == 0)
		{
			maxx = c[b[i].start].x;
		}
		if (c[b[i].start].x == c[b[i].end].x)
		{
			if (c[b[i].start].x > maxx)
			{
				maxx = c[b[i].start].x;
			}
			all.push_back(b[i]);
		}
	}
	for (int i = 0; i < all.size(); i++)
	{
		if (c[all[i].start].x == maxx)
		{
			rightbounds.push_back(all[i]);
		}
	}
	return rightbounds;

}
vector<edge> getupperbounds(vector<edge> b, vector<vertex> c)
{
	vector<edge> lowerbounds;
	vector<edge> upperbounds;
	vector<edge> all;
	double maxx;
	for (int i = 0; i < b.size(); i++)
	{
		if (i == 0)
		{
			maxx = c[b[i].start].y;
		}
		if (c[b[i].start].y == c[b[i].end].y)
		{
			if (c[b[i].start].x > maxx)
			{
				maxx = c[b[i].start].y;
			}
			all.push_back(b[i]);
		}
	}
	for (int i = 0; i < all.size(); i++)
	{
		if (c[all[i].start].y== maxx)
		{
			upperbounds.push_back(all[i]);
		}
	}
	return upperbounds;

}
vector<edge> getlowerbounds(vector<edge> b, vector<vertex> c)
{
	vector<edge> lowerbounds;
	vector<edge> upperbounds;
	vector<edge> all;
	double minx;
	for (int i = 0; i < b.size(); i++)
	{
		if (i == 0)
		{
			minx = c[b[i].start].y;
		}
		if (c[b[i].start].y == c[b[i].end].y)
		{
			if (c[b[i].start].x < minx)
			{
				minx = c[b[i].start].y;
			}
			all.push_back(b[i]);
		}
	}
	for (int i = 0; i < all.size(); i++)
	{
		if (c[all[i].start].y == minx)
		{
			lowerbounds.push_back(all[i]);
		}
	}
	return lowerbounds;

}
vector<vector<float>> getni(vector<edge> b, vector<vertex> c)
{
	float nx, ny;
	float l;
	vector < vector<float> > ni;
	vector<float> ins;
	ins.resize(3);
	for (int i = 0; i < b.size(); i++)
	{
		nx = c[b[i].end].x - c[b[i].start].x;
		ny= c[b[i].start].y - c[b[i].end].y;
		l = sqrt(pow(nx, 2) + pow(ny, 2));
		nx = nx / l;
		ny = ny / l;
		int index;
		for (int j = 0; j < c[b[i].start].connections.size(); j++)
		{
			bool check=true;
			for (int k = 0; k < c[b[i].end].connections.size(); k++)
			{
				if (c[b[i].start].connections[j] == c[b[i].end].connections[k])
				{
					index = c[b[i].start].connections[j];
					check = false;
					break;

				}
			}
			if (!check)
			{
				break;
			}
		}
		float testx =  c[b[i].end].x- c[index].x;
		float testy= c[b[i].end].y - c[index].y;
		if (testx * nx + testy * ny < 0)
		{
			nx = -nx;
			ny = -ny;

		}
		ins[0] = nx;
		ins[1] = ny;
		ins[2] = l;
		ni.push_back(ins);


	}
	return ni;
}
Eigen::SparseMatrix<float> getp(vector<edge> b, vector<vertex> c,float P)
{
	Eigen::SparseMatrix<float> f(2*c.size(),1);
	vector<Eigen::Triplet<float>> triplets;
	vector<vector<float>> ni = getni(b, c);
	for (int i = 0; i < b.size(); i++)
	{
		float l = ni[i][2];
		/**global_F[2 * data[i][1], 0] += (P * n_vec[i][0] * l)
			global_F[2 * data[i][1] + 1, 0] += (P * n_vec[i][1] * l)
			global_F[2 * data[i][2], 0] += (P * n_vec[i][0] * l)
			global_F[2 * data[i][2] + 1, 0] += (P * n_vec[i][1] * l)*/
		triplets.push_back(Eigen::Triplet<float>(2 * b[i].start,0, P * ni[i][0] * l));
		triplets.push_back(Eigen::Triplet<float>(2 * b[i].start+1, 0, P * ni[i][1] * l));
		triplets.push_back(Eigen::Triplet<float>(2 * b[i].end, 0, P * ni[i][0] * l));
		triplets.push_back(Eigen::Triplet<float>(2 * b[i].end+1, 0, P * ni[i][1] * l));
		


	}
	f.setFromTriplets(triplets.begin(), triplets.end());
	f.makeCompressed();
	return f;


}
Eigen::SparseMatrix<float> generatek(Eigen::Matrix3f D,vector<element> el, vector<vertex> coord)
{
	Eigen::SparseMatrix<float> k(2*coord.size(), 2 * coord.size());
	vector<Eigen::Triplet<float>> triplets;
	//vector<Eigen::Triplet<float>>& trp =triplets;
	for (std::vector<element>::iterator it = el.begin(); it != el.end(); ++it)
	{
		triplets=it->CalculateStiffnessMatrix(D, triplets, coord);
	}
	/*for (int i = 0; i < triplets.size(); i++)
		cout << triplets[i].row() << " " << triplets[i].col() << " " << triplets[i].value()<<endl;*/
	k.setFromTriplets(triplets.begin(), triplets.end());
	return k;
}
void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index)
{
	if (it.row() == index || it.col() == index)
	{
		it.valueRef() = it.row() == it.col() ? 1.0f : 0.0f;
		
	}
}

void applyconstraints(vector<edge> b, Eigen::SparseMatrix<float> &K, Eigen::SparseMatrix<float> &F)
{
	//x-val->0;
	//y-va->+1
	for (int k = 0; k < K.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<float>::InnerIterator it(K, k); it; ++it)
		{
			for (int i = 0; i < b.size(); i++)
			{
				
				SetConstraints(it, 2 * b[i].start); 
				SetConstraints(it, 2 * b[i].start+1);
				SetConstraints(it, 2 * b[i].end);
				SetConstraints(it, 2 * b[i].end+1);
			}
		}
	}
	for (int k = 0; k < F.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<float>::InnerIterator it(F, k); it; ++it)
		{
			for (int i = 0; i < b.size(); i++)
			{
				if (2 * b[i].start == k|| 2 * b[i].start + 1==k|| 2 * b[i].end==k|| 2 * b[i].end + 1==k)
				{
					it.valueRef() = 0;
				}
			}
		}
	}

	
}
int main()
{
	string txt;
	//fptr = fopen("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug.k");
	string filename = "C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\solution python\\Kirsch_mesh_hirarch_size_small.k";
	//ifstream newfile("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug.k");
	//string filename = "C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug.k";
	ifstream newfile(filename);
	ofstream nf("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug_arr.txt");
	//newfile.open("C:\\Users\\maria\\Documents\\univ\\2022\\fall\\prak\\kiesh_treug.k")
	vector<vertex> arr_Coord;
	vector<element> arr_elem;
	vector<edge>arr_edge;
	std::vector< Eigen::Triplet<int> > triplets;
	//матрица инцендентности
	while (std::getline(newfile, txt) )
	{
		if (txt == "*NODE")
		{
			break;
		}
	}
	std::getline(newfile, txt);
	
	while (std::getline(newfile, txt) )
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
	std::getline(newfile, txt);
	std::getline(newfile, txt);
	std::getline(newfile, txt);
	
	
	while (std::getline(newfile, txt))
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
			arr_Coord[n].in_element_id.resize(arr_Coord[n].in_element_id.size () + 1);
			arr_Coord[n].in_element_id[arr_Coord[n].in_element_id.size() - 1] = arr_elem.size()-1;
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
				arr_Coord[n].in_element_id.resize(arr_Coord[n].in_element_id.size() + 1);
				arr_Coord[n].in_element_id[arr_Coord[n].in_element_id.size() - 1] = arr_elem.size() - 1;
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
						//cout << j;
							//j=
					}
				//vector<int>::iterator it=find(arr_Coord[i].connections.begin(), arr_Coord[i].connections.end(), j);
				

				//if (it == arr_Coord[i].connections.end())
				//{
					triplets.push_back(Eigen::Triplet<int>(i, j, 1));
					triplets.push_back(Eigen::Triplet<int>(j, i, 1));
					
					arr_Coord[i].connections.resize(arr_Coord[i].connections.size() + 1);
					arr_Coord[i].connections[arr_Coord[i].connections.size() - 1] = j;
					arr_Coord[j].connections.resize(arr_Coord[j].connections.size() + 1);
					arr_Coord[j].connections[arr_Coord[j].connections.size() - 1] = i;

				//}
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
		if (mindegree > arr_Coord[i].degree && i!=0)
		{
			mindegree = arr_Coord[i].degree;
			idmind = i;
		}
		
		/*for (int k = 0; k < arr_Coord[i].connections.size(); k++)
		{
			triplets.push_back(Eigen::Triplet<int>(i, arr_Coord[i].connections[k], 1));
		}*/
	}
	adjacency.setFromTriplets(triplets.begin(), triplets.end());
	//adjacency.insert(0, 0);
	adjacency.makeCompressed();
	for (int k = 0; k < adjacency.outerSize(); ++k)
		for (Eigen::SparseMatrix<int>::InnerIterator it(adjacency, k); it; ++it)
		{
			if (it.row() <= it.col())
			{
				arr_edge.resize(arr_edge.size() + 1);
				arr_edge[arr_edge.size() - 1].in_element = it.value();
				arr_edge[arr_edge.size() - 1].start = it.row();
				arr_edge[arr_edge.size() - 1].end = it.col();
			}
			//cout << it.value() << ' ';
			/*cout << it.row() << ' ';   // row index
			cout<<it.col() << ' ';  // col index (here it is equal to k)
			cout << endl;*/
			//cout<<it.index() << ' '; // inner index, here it is equal to it.row()
		}
	std::cout << "Start" << endl;
	//draw_matrix(adjacency,"before");
	//draw_coord(arr_Coord);
	draw_edges(arr_edge, arr_Coord,"setka");
	//Cuthill_Mckee(idmind,adjacency,arr_Coord,arr_elem,arr_edge);
	adjacency = build_matrix(arr_Coord);
	//draw_matrix(adjacency,"after");
	std::cout << "Nodes: " << arr_Coord.size() << endl << "Elememts " << arr_elem.size();
	std::cout << std::endl;//<< adjacency;
	vector<edge> bounds;
	bounds = getbounds(arr_edge);
	//draw_edges(bounds, arr_Coord);
	
	vector<edge> leftbounds;
	vector<edge> rightbounds;
	vector<edge> upbounds;
	vector<edge> downbounds;
	leftbounds = getleftbounds(bounds,arr_Coord);
	rightbounds = getrightbounds(bounds, arr_Coord);
	upbounds = getupperbounds(bounds, arr_Coord);//need to change their palaces
	downbounds = getlowerbounds(bounds, arr_Coord);
	//draw_edges(downbounds, arr_Coord);

	

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

	Eigen::SparseMatrix<float> K;
	Eigen::SparseMatrix<float> F;
	K = generatek(D, arr_elem, arr_Coord);
	F = getp(upbounds, arr_Coord, 100000000);
	applyconstraints(leftbounds, K, F);
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver(K);
	Eigen::VectorXf displacements = solver.solve(F);
	//cout << displacements;
	//sigma=dbu
	//y=0;
	double xmin=100;
	double xmax = 0;

	for (int i = 0; i < downbounds.size(); i++)
	{
		if (xmax < arr_Coord[downbounds[i].start].x)
		{
			xmax = arr_Coord[downbounds[i].start].x;
		}

		if (xmax < arr_Coord[downbounds[i].end].x)
		{
			xmax = arr_Coord[downbounds[i].end].x;
		}
		if (xmin > arr_Coord[downbounds[i].start].x)
		{
			xmin = arr_Coord[downbounds[i].start].x;
		}
		if (xmin > arr_Coord[downbounds[i].end].x)
		{
			xmin = arr_Coord[downbounds[i].end].x;
		}
	}
	int index=closest(xmin, 0, arr_Coord);
	cout << index << " " << arr_Coord[index].x << " " << arr_Coord[index].y << endl;
	int colst = 1000;//stepamount
	double h = (xmax - xmin) / colst;
	//double x = xmin;
	vector<float> x;
	for (int i = 0; i < colst; i++)
	{
		x.push_back(xmin + i * h);
	}
	
		vector<float> y;
	y=analytical(xmin, 0, xmax, 0, h);
	draw_graph(x, y);




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