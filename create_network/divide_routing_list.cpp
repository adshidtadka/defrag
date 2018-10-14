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

	int path_read[5];
	ifs_mixed >> path_read[0] >> path_read[1] >> path_read[2] >> path_read[3] >> path_read[4];
	int s_prev = 0;
	int path_num = 0;
	while(path_read[0] >= s_prev){
		s_prev = path_read[0];
		for (int i = 0; i < 4; ++i)
		{
			path[path_num][i] = path_read[i];
		}
		ifs_mixed >> path_read[0] >> path_read[1] >> path_read[2] >> path_read[3] >> path_read[4];
		path_num++;
	}
	path_num--;
	ifs_mixed.close();

	return path_num;
}

int copyToPathSameSD(int path[10000][4], int path_same_sd[20][4], int path_count, int path_same_sd_count){
	for (int j = 0; j < 4; ++j)
	{
		path_same_sd[path_same_sd_count][j] = path[path_count][j];
	}
	path_same_sd_count++;
	return path_same_sd_count;
}

int writeRoot(int path_same_sd[20][4], int path_same_sd_count, int first_dest_prim, ofstream &ofs)
{
	int start_node = path_same_sd[0][0];		//source of root
	int end_node   = path_same_sd[0][1];		//destination of root
	int next_source;							//next source of link

	for (int i = 0; i < path_same_sd_count; ++i)
	{
		if (path_same_sd[i][2] == start_node && path_same_sd[i][3] != first_dest_prim)
		{	
			// save next source node
			next_source = path_same_sd[i][3];
			first_dest_prim = path_same_sd[i][3];

			// write first prim link
			for (int j = 0; j < 4; ++j)
			{
				ofs << path_same_sd[i][j] << " ";
			}

			// delete path
			path_same_sd[i][2] = -1;

			ofs << " 1" << endl;
	
	
			//write the other links
			for (int j = 0; j < path_same_sd_count; ++j)
			{
				if (next_source == end_node)
				{
					// finish write the other links
					break;
				}
				if (path_same_sd[j][2] == next_source)
				{
					// save next source node
					next_source = path_same_sd[j][3];

					// write a link
					for (int k = 0; k < 4; ++k)
					{
						ofs << path_same_sd[j][k] << " ";

					}
					ofs << " 1" << endl;

					// delete path
					path_same_sd[j][2] = -1;
	
					//reset search index
					j = -1;
				}
			}
			// finish write the other links
			break;
		}
	}
	return first_dest_prim;
}

int writeRoots(int path[10000][4], int argc, char* argv[0], int path_num){

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
	// write first line
	ofs_prim << "f :=" << endl;
	ofs_back << "f :=" << endl;

	int path_same_sd[20][4];
	int path_same_sd_count = 0;
	path_same_sd[0][0]  = -1;
	for (int i = 0; i <= path_num; ++i)
	{
		// next candidate of s-d
		if (path_same_sd[0][0]  == -1)
		{
			path_same_sd_count = copyToPathSameSD(path, path_same_sd, i, path_same_sd_count);
			continue;
		}

		// same s-d
		if (path[i][0] == path_same_sd[path_same_sd_count -1][0] && path[i][1] == path_same_sd[path_same_sd_count -1][1])
		{
			path_same_sd_count = copyToPathSameSD(path, path_same_sd, i, path_same_sd_count);
			continue;
		} else {

			int first_dest_prim = -1;

			// write path_prim
			first_dest_prim = writeRoot(path_same_sd, path_same_sd_count, first_dest_prim, ofs_prim);

			// write path_back
			writeRoot(path_same_sd, path_same_sd_count, first_dest_prim, ofs_back);

			//initialize
			path_same_sd[0][0]  = -1;
			path_same_sd_count = 0;
			i--;
			continue;
		}
	}
	ofs_prim << ";" << endl;
	ofs_back << ";" << endl;
	ofs_prim.close();
	ofs_back.close();


	return 0;
}

int main(int argc, char* argv[])
{
	int mixed_path[10000][4];
	int path_num;

	path_num = readMixedList(mixed_path, argc, &argv[0]);
	writeRoots(mixed_path, argc, &argv[0], path_num);
	 
	return 0;
}