#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

int readMixedList(int path[10000][4], int argc, char* argv[0]){

	ifstream ifs_mixed;

	ifs_mixed.open (argv[1]);
	if (!ifs_mixed){
		cout <<"Cannot open mixed list" << endl;
		return 1;
	}

	ifs_mixed.ignore(INT_MAX,'=');

	int path_read[5] = {0, 0, 0, 0, 0};
	int s_prev = 0;
	int path_num = 0;
	while(path_read[0] >= s_prev){
		s_prev = path_read[0];
		for (int i = 0; i < 4; ++i)
		{
			path[path_num][i] = path_read[i];
		}
		ifs_mixed >> path_read[0] >> path_read[1] >> path_read[2] >> path_read[3] >> path_read[4];
		// cout << "path_read[0] = " << path_read[0] << ", path_read[1] = " << path_read[1] << ", path_read[2] = " << path_read[2] << ", path_read[3] = " << path_read[3] << ", path_read[4] = " << path_read[4] <<  ", s_prev = " << s_prev <<  endl;
		path_num++;
	}
	path_num--;
	cout << "path_num = " << path_num << endl;
	ifs_mixed.close();

	return path_num;
}

int insertPath(int path[10000][4], int path_same_sd[20][4], int key_path, int key_same){
	for (int j = 0; j < 4; ++j)
	{
		path_same_sd[key_same][j] = path[key_path][j];
	}
	key_same++;
	return key_same;
}

// int addDividedList(int path[10000][4], int path_same_sd[20][4], int key){
// 	for (int i = 0; i < key; ++i)
// 	{
		
// 	}
// 	return 0;
// }

int writeDevidedList(int path[10000][4], int argc, char* argv[0], int path_num){

	ofstream ofs_prim;
	ofstream ofs_back;

	ofs_prim.open (argv[2]);
	ofs_back.open (argv[3]);
	if (!ofs_prim){
		cout <<"Cannot open prim list" << endl;
		return 1;
	}

	if (!ofs_back){
		cout <<"Cannot open back list" << endl;
		return 1;
	}

	ofs_prim << "f =:" << endl;
	ofs_back << "f =:" << endl;

	//devide two mixed list
	int path_same_sd[20][4];
	int key = 0;
	path_same_sd[0][0]  = -1;
	for (int i = 0; i < path_num; ++i)
	{
		// next candidate of s-d
		if (path_same_sd[0][0]  == -1)
		{
			key = insertPath(path, path_same_sd, i, key);
			continue;
		}

		// same s-d
		if (path[i][0] == path_same_sd[key -1][0] && path[i][1] == path_same_sd[key -1][1])
		{
			key = insertPath(path, path_same_sd, i, key);
			continue;
		} else {
			// add path_same_sd
			// addDividedList(path, path_same_sd, key);

			//initialize
			int path_same_sd[20][4] = {{}};
			path_same_sd[0][0]  = -1;
			key = 0;
			i--;
			continue;
		}
	}
	 ofs_prim.close();
	 ofs_back.close();


	return 0;
}

int main(int argc, char* argv[])
{
	int mixed_path[10000][4];
	int path_num;

	path_num = readMixedList(mixed_path, argc, &argv[0]);

	// writeDevidedList(mixed_path, argc, &argv[0], path_num);


	 
	return 0;
}