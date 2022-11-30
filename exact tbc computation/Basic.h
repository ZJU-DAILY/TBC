#pragma once
#include <string>
#include <utility>
#include <vector>
using namespace std;
template <class T>
struct Node {
    int seq;
    string identity;
    vector<T> adj_edge;
    Node() = default;
    Node(int _seq, string _id,  vector<T> _adj) :
        seq(_seq), identity(move(_id)),  adj_edge(move(_adj)) {}
};
struct Edge {
    int source_seq;
    int destination_seq;
    int timestamp;
    Edge() = default;
    Edge(int _src, int _des, int _timestamp) :
        source_seq(_src), destination_seq(_des), timestamp(_timestamp) {}
};
struct DataNode : Node<Edge> {
    DataNode() = default;
    DataNode(int _seq, string _id,  vector<Edge> _adj) :
        Node(_seq, move(_id), move(_adj)) {}
};
struct Tran_Edge {
    int source_seq;
    int destination_seq;

    Tran_Edge() = default;
    Tran_Edge(int _src, int _des) :
        source_seq(_src),  destination_seq(_des) {}
};
struct Tran_Node
{
    int seq;
    int timestamp;
    int identity;
    vector<Tran_Edge> adj_edge;
    Tran_Node() = default;
    Tran_Node(int _seq, int _timestamp, int _identity, vector<Tran_Edge> _adj) :
        seq(_seq), timestamp(_timestamp), identity(_identity), adj_edge(_adj) {}

};
struct Cmp_Edge {
    int source_seq;
    int destination_seq;
    Cmp_Edge() = default;
    Cmp_Edge(int _src,  int _des) :
        source_seq(_src),  destination_seq(_des) {}
};
struct Cmp_Node
{
    int seq;
    vector<int> timestamps;
    int identity;
    vector<Cmp_Edge> adj_edge;
    Cmp_Node() = default;
    Cmp_Node(int _seq, vector<int> _timestamps, int _identity, vector<Cmp_Edge> _adj) :
        seq(_seq), timestamps(_timestamps), identity(_identity), adj_edge(_adj) {}

};

//struct ParentScatterPoint {
//    pair<int, int> parent;
//    int num;
//};