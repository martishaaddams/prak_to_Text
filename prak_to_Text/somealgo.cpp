#include"somealgo.h"
void renumerate(std::vector<vertex>& coord, std::map<int, int> renum)
{
	for (int i = 0; i < coord.size(); i++)
	{
		coord[i].renum(renum);
	}
	std::vector<vertex> tmp;
	tmp.resize(coord.size());
	tmp = coord;
	//vertex temp;
	for (int i = 0; i < coord.size(); i++)
	{
		//temp = coord[i];
		coord[i] = tmp[renum[i]];
		//coord[i] = temp;
	}
	
}
void renumerate(std::vector<vertex>& coord, std::vector<int> renum)
{
	for (int i = 0; i < coord.size(); i++)
	{
		coord[i].renum(renum);
	}
	std::vector<vertex> tmp;
	tmp.resize(coord.size());
	tmp = coord;
	for (int i = 0; i < coord.size(); i++)
	{
		tmp[renum[i]] = coord[i];
	}
	coord = tmp;
}
void renumerate(std::vector<element>& el, std::map<int, int> renum)
{
	for (int i = 0; i < el.size(); i++)
	{
		el[i].renum(renum);
	}
}
void renumerate(std::vector<edge>& el, std::vector<int> renum)
{
	for (int i = 0; i < el.size(); i++)
	{
		el[i].renumerate(renum);
	}
}
void print(std::vector<int> a)
{
	for (auto i = 0; i < a.size(); i++)
	{
		std::cout << a[i] << ' ';
	}
	std::cout << std::endl;
}
bool compare_vertex_deg(int v1, int v2, std::vector<vertex> coord)
{
	return(coord[v1].degree < coord[v2].degree);
	//return(v1.x < v2.x);
}
void swap(int* a, int* b)
{
	int t = *a;
	*a = *b;
	*b = t;
}
int partition(std::vector<int>& arr, int low, int high, std::vector<vertex> coord)
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
void quickSort(std::vector<int>& arr, int low, int high, std::vector<vertex> coord)
{
	if (low < high) {
		/* pi is partitioning index, arr[p] is now
		at right place */
		int pi = partition(arr, low, high, coord);

		// Separately sort elements before
		// partition and after partition
		quickSort(arr, low, pi - 1, coord);
		quickSort(arr, pi + 1, high, coord);
	}
}
int findIndex(std::vector<std::pair<int, int> > a, int x)
{
	for (int i = 0; i < a.size(); i++)
		if (a[i].first == x)
			return i;
	return -1;
}
std::map<int, int> Cuthill_Mckee(int indexofvert, Eigen::SparseMatrix<int> mat, std::vector<vertex>& coord, std::vector<element>& el,std::vector<edge> &ed)
{
	Eigen::SparseMatrix<int> res(mat.rows(), mat.cols());

	std::map<int, int> renum;
	std::queue<int> qu;
	std::vector<int> rn;
	rn.resize(coord.size());
	//std::vector<std::pair<int,int>> notVisited;

	//for (int i = 0; i < coord.size(); i++)
	//{
	//	notVisited.push_back(std::make_pair(i, coord[i].degree));
	//}
	//while (notVisited.size())
	//{
	//	int minNodeIndex = 0;
	//	//this is useless if graph connected
	//	for (int i = 1; i < notVisited.size(); i++)
	//		if (notVisited[i].second < notVisited[minNodeIndex].second)
	//			minNodeIndex = i;
	//	qu.push(notVisited[minNodeIndex].first);

	//	notVisited.erase(notVisited.begin() + findIndex(notVisited, notVisited[qu.front()].first));
	//	int k = 1;
	//	while (!qu.empty())
	//	{
	//		if (renum.find(qu.front()) == renum.end())
	//		{
	//			std::vector<int> toSort;;
	//			for (int i = 0; i < coord[qu.front()].connections.size(); i++)
	//			{
	//				if (findIndex(notVisited, i) != -1)
	//				{
	//					toSort.push_back(i);
	//					notVisited.erase(notVisited.begin() + findIndex(notVisited, i));
	//				}
	//			}
	//			quickSort(toSort, 0, toSort.size() - 1, coord);
	//			for (int i = 0; i < toSort.size(); i++)
	//				qu.push(toSort[i]);
	//			renum[qu.front()] = k;
	//			k++;

	//		}
	//		qu.pop();
	//	}
	//}
	////renum[indexofvert] = 1;
	qu.push(indexofvert);

	int k = 1;
	while (!qu.empty())
	{
		if (renum.find(qu.front()) == renum.end())
		{
			
			quickSort(coord[qu.front()].connections, 0, coord[qu.front()].connections.size() - 1, coord);
			
			//print(coord[qu.front()].connections);
			/*for (int i = 0; i < coord[qu.front()].connections.size(); i++)
			{
				std::cout << coord[coord[qu.front()].connections[i]].degree << ' ';
			}*/
			//std::cout << std::endl;
			for (auto i = 0; i < coord[qu.front()].connections.size(); i++)
			{
				qu.push(coord[qu.front()].connections[i]);
			}

			renum[qu.front()] = k;
			rn[qu.front()] =k;
			//renumx[k] = qu.front();
			k++;
			//std::cout << qu.front() << std::endl;
		}
		qu.pop();
	}
	//for (auto it = 0; it <rn.size(); it++)
	//{
	//	//std::cout << it->first << '->' << it->second << std::endl;
	//	std::cout << it << ' '<<rn[it]<<' '<<std::endl;

	//}
	renumerate(el, renum);
	renumerate(coord, rn);
	renumerate(ed, rn);
	return renum;
}
