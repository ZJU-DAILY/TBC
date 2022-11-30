#pragma once
#include "Basic.h"
#include <map>
using namespace std;
//time instance graph is compressed
struct CompressGraph 
{
	vector<Tran_Node> tran_data_nodes;
	vector<Cmp_Node> cmp_data_nodes;
	map<int, int> mapping;

	CompressGraph() = default;
	CompressGraph(vector<Tran_Node> _tran_data_nodes) :
		tran_data_nodes(_tran_data_nodes) {}
};