#pragma once
#include <unordered_map>
#include <map>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <algorithm>
#include <functional>
#include <queue>
#include <numeric>
#include <unordered_set>
#include <string>
#include "Basic.h"
using namespace std;
namespace std
{
    template<>
    struct hash<pair<pair<int, int>, int>> {
    public:
        size_t operator()(const pair<pair<int, int>, int>& p)const {
            return hash<int>()(p.first.first) ^ hash<int>()(p.first.second) ^ hash<int>()(p.second);
        }
    };
    template<>
    struct hash<pair<int, int>> {
    public:
        size_t operator()(const pair<int, int>& p)const {
            return hash<int>()(p.first) ^ hash<int>()(p.second);
        }
    };
}
struct DataGraph {
    vector<DataNode> graph_data;
    unordered_map<string, int> identity_to_seq;
    unordered_map<int, string> seq_to_identity;
    map<int, int> timestamps;
    unordered_set<pair<pair<int, int>, int>> edge_set;
    int max_ts = 0;
    void loader(string file_name) {
        int n = 0;
        int m = 0;
        graph_data.clear();
        identity_to_seq.clear();
        ifstream infile(file_name.c_str());
        if (!infile) {
            cout << "open data-graph " << file_name << " filed!\n";
            return;
        }
        string line;
        while (getline(infile, line)) 
        {
            if (line.empty()) break;
            istringstream instr(line);
            vector<string> id(2);
            string timestamp;
            instr >> id[0] >> id[1] >> timestamp;

            string tmp = id[0] + ' ' + id[1] + ' ' + timestamp;
            m++;
            for (int i = 0; i < 2; ++i) {
                if (identity_to_seq.find(id[i]) == identity_to_seq.end()) {
                    n++;
                    DataNode new_node(graph_data.size(), id[i], {});
                    graph_data.push_back(new_node);
                    identity_to_seq[new_node.identity] = new_node.seq;
                    seq_to_identity[new_node.seq] = new_node.identity;
                }
            }
            if (id[0] == id[1]) continue;
            if (!edge_set.count({ {identity_to_seq[id[0]], identity_to_seq[id[1]]},stoi(timestamp) }))
            {
                edge_set.insert({ {identity_to_seq[id[0]], identity_to_seq[id[1]]},stoi(timestamp) });
                Edge new_edge(identity_to_seq[id[0]], identity_to_seq[id[1]], stoi(timestamp));
                graph_data[identity_to_seq[id[0]]].adj_edge.push_back(new_edge);
                timestamps.insert({ new_edge.timestamp ,0 });
            }
        }
        //开始时间和结束时间升序排序
        for (auto& node : graph_data) {
            sort(node.adj_edge.begin(), node.adj_edge.end(), [](const Edge& lhs, const Edge& rhs) {
                return lhs.timestamp < rhs.timestamp;
                });
        }
        //        cout << graph_data.size() << '\n';
        cmp_ts();
        cout << "set data-graph " << file_name << " successfully!\n";
        cout << "----project name:Tranform-----\n";
        cout << "node size:" << n << endl;
        cout << "edge size:" << m << endl;
        //cout << "rep_edges: " << rep_edges << endl;
        infile.close();

    }
    void cmp_ts()
    {
        int ts = 1;
        for (auto& timestamp_pair : timestamps)
        {
            timestamp_pair.second = ts;
            ts++;
        }
        max_ts = ts - 1;
        cout << "maximal timestamp: " << max_ts << endl;
        for (int i = 0; i < graph_data.size(); i++)
        {
            for (int j = 0; j < graph_data[i].adj_edge.size(); j++)
            {
                Edge& e = graph_data[i].adj_edge[j];
                e.timestamp = timestamps[e.timestamp];
            }
        }
    }
};