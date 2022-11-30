#pragma once
#include "TransformGraph.h"
#include "CompressGraph.h"
#include "DataGraph.h"
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <omp.h>
#include <climits>
struct Transformer
{
	DataGraph data_graph;
	TransformGraph tran_graph;
	vector<DataNode>& data_nodes = data_graph.graph_data;
	vector<Tran_Node>& tran_nodes = tran_graph.Tran_Data_Nodes;
	unordered_map<int, string>& seq_to_identity = data_graph.seq_to_identity;
	unordered_map<string, int>& identity_to_seq = data_graph.identity_to_seq;
	CompressGraph cmp_graph;
	unordered_map<int, int> tran_root_seqs;
	int cmp_edge_size = 0;
	int tran_edge_size = 0;
	unordered_set<int> needcreate;
	Transformer(string _filename)
	{
		tran_graph = TransformGraph();
		data_graph = DataGraph();
		data_graph.loader(_filename);
	}
	//orgin graph -> time instance graph
	void doTransform(bool isStrict)
	{
		int strict = isStrict ? 0 : 1;
		unordered_set<pair<int, int>>existed_sp;
		unordered_map<int, vector<int>> tran_edges;
		unordered_map<pair<int, int>, int> mapping;
		for (int i = 0; i < data_nodes.size(); i++)
		{
			unordered_set<int> tmp_adj_sp;
			Tran_Node tran_node = Tran_Node(tran_nodes.size(), 0, data_nodes[i].seq, {});
			tran_root_seqs[i] = tran_node.seq;
			mapping[{i, 0}] = tran_node.seq;
			vector<pair<int, int>> tmp_edges;
			for (int j = 0; j < data_nodes[i].adj_edge.size(); j++)
			{
				pair<int, int> sp = { data_nodes[i].adj_edge[j].destination_seq, data_nodes[i].adj_edge[j].timestamp };
				if (existed_sp.count(sp))
				{
					int tran_sp = mapping[sp];
					tmp_adj_sp.insert(tran_sp);
				}
				else tmp_edges.push_back({ i,j });
			}
			for (auto& tmp_tran_sp : tmp_adj_sp)
			{
				tran_node.adj_edge.push_back(Tran_Edge(tran_node.seq, tmp_tran_sp));
				tran_edge_size++;
			}
			tran_nodes.push_back(tran_node);

			if (tmp_edges.size() > 0)
			{
				queue<pair<int, int>> q;
				for (auto& e : tmp_edges)
				{
					int des_seq = data_nodes[e.first].adj_edge[e.second].destination_seq;
					int ts = data_nodes[e.first].adj_edge[e.second].timestamp;
					pair<int, int> sp_ = { des_seq,ts };
					Tran_Node node = Tran_Node(tran_nodes.size(), sp_.second, data_nodes[des_seq].seq, {});
					tran_nodes.push_back(node);
					mapping[{des_seq, ts}] = node.seq;
					tran_edges[tran_node.seq].push_back(node.seq);
					q.push({ des_seq,node.timestamp });
					existed_sp.insert(sp_);
				}
				while (!q.empty())
				{
					pair<int, int> sp = q.front();
					q.pop();
					for (int j = 0; j < data_nodes[sp.first].adj_edge.size(); j++)
					{
						pair<int, int> e = { sp.first,j };
						pair<int, int> sp_ = { data_nodes[sp.first].adj_edge[j].destination_seq,data_nodes[sp.first].adj_edge[j].timestamp };
						if (sp.second < sp_.second + strict)
						{
							if (existed_sp.count(sp_))
							{
								int tran_sp = mapping[sp];
								int tran_sp_ = mapping[sp_];
								tran_edges[tran_sp].push_back(tran_sp_);
							}
							else
							{
								Tran_Node node = Tran_Node(tran_nodes.size(), sp_.second, data_nodes[sp_.first].seq, {});
								tran_nodes.push_back(node);
								mapping[sp_] = node.seq;
								tran_edges[mapping[sp]].push_back(node.seq);
								q.push(sp_);
								existed_sp.insert(sp_);
							}
						}
					}
				}
			}
		}

		for (auto edge_info : tran_edges)
		{
			int u = edge_info.first;
			for (int i = 0; i < edge_info.second.size(); i++)
			{
				int u_ = edge_info.second[i];
				Tran_Node& node = tran_nodes[u];
				Tran_Node& node_ = tran_nodes[u_];
				node.adj_edge.push_back(Tran_Edge(node.seq, node_.seq));
				tran_edge_size++;
			}
		}
		int tran_node_size = tran_nodes.size();

		unordered_map<int, int> inverse_edge_min_ts;
		for (int i = 0; i < data_nodes.size();i++)
		{
			for (int j = 0; j < data_nodes[i].adj_edge.size(); j++) 
			{
				int to = data_nodes[i].adj_edge[j].destination_seq;
				if (inverse_edge_min_ts.count(to))
				{
					inverse_edge_min_ts[to] = min(inverse_edge_min_ts[to], data_nodes[i].adj_edge[j].timestamp);
				}
				else
				{
					inverse_edge_min_ts[to] = data_nodes[i].adj_edge[j].timestamp;
				}
			}
		}

		for (int i = 0; i < data_nodes.size(); i++) 
		{
			for (int j = 0; j < data_nodes[i].adj_edge.size(); j++) 
			{
				int timestamp = data_nodes[i].adj_edge[j].timestamp;
				if (!inverse_edge_min_ts.count(i)) 
				{
					needcreate.insert(i);
					break;
				}
				if (inverse_edge_min_ts.count(i) && inverse_edge_min_ts[i] > timestamp) 
				{
					needcreate.insert(i);
					break;
				}
			}
		}
		for (auto& p : tran_root_seqs) 
		{
			int i = p.first;
			int tran_seq = p.second;
			if (!needcreate.count(i)) 
			{

				tran_node_size--;
				tran_edge_size = tran_edge_size - tran_nodes[tran_seq].adj_edge.size();
			}
		}

	}
	// vertex and edge compression
	void doCompress()
	{
		vector<Cmp_Node>& cmp_nodes = cmp_graph.cmp_data_nodes;
		unordered_map<int, unordered_map<int, vector<int>>> in_edges;
		unordered_map<int, unordered_map<int, vector<int>>> out_edges;
		unordered_set<int> exist_tran_seq;
		unordered_map<int, vector<int>> same_id_diff_seq;
		unordered_map<int, int> tran_seq_to_cmp_seq;
		unordered_set<int> vertical_tran_seq;
		for (int j = 0; j < tran_nodes.size(); j++)
		{
			int tran_seq = j;
			if (tran_nodes[tran_seq].timestamp != 0)
				same_id_diff_seq[tran_nodes[tran_seq].identity].push_back(tran_seq);
			for (int k = 0; k < tran_nodes[tran_seq].adj_edge.size(); k++)
			{
				Tran_Edge& e = tran_nodes[tran_seq].adj_edge[k];
				out_edges[tran_nodes[e.source_seq].identity][e.source_seq].push_back(e.destination_seq);
				in_edges[tran_nodes[e.destination_seq].identity][e.destination_seq].push_back(e.source_seq);
			}
		}
		for (int j = 0; j < tran_root_seqs.size(); j++)
		{
			int tran_seq = tran_root_seqs[j];
			Cmp_Node cmp_node = Cmp_Node(cmp_nodes.size(), { 0 }, tran_nodes[tran_seq].identity, {});
			cmp_nodes.push_back(cmp_node);
			tran_seq_to_cmp_seq.insert({ tran_seq,cmp_node.seq });
		}
		unordered_set<int> tran_seq_exist;
		for (auto& same_id_info : same_id_diff_seq)
		{
			if (same_id_info.second.size() == 1)
			{
				int tran_seq = same_id_info.second.front();
				Cmp_Node cmp_node = Cmp_Node(cmp_nodes.size(), { tran_nodes[tran_seq].timestamp }, tran_nodes[tran_seq].identity, {});
				cmp_nodes.push_back(cmp_node);
				tran_seq_to_cmp_seq.insert(make_pair(tran_seq, cmp_node.seq));
			}
			else
			{
				int next_n = 0;
				for (int n = 0; n < same_id_info.second.size();)
				{
					for (int m = n + 1; m < same_id_info.second.size(); m++)
					{
						int pair1 = same_id_info.second[n];
						int pair2 = same_id_info.second[m];
						vector<int>& in_edge_n = in_edges[tran_nodes[pair1].identity][pair1];
						vector<int>& out_edge_n = out_edges[tran_nodes[pair1].identity][pair1];
						vector<int>& in_edge_m = in_edges[tran_nodes[pair2].identity][pair2];
						vector<int>& out_edge_m = out_edges[tran_nodes[pair2].identity][pair2];

						bool in = false;
						bool out = false;
						
						if (in_edge_n.size() == in_edge_m.size())
						{
							int x = 0;
							for (x; x < in_edge_n.size(); x++)
							{
								if (in_edge_n[x] != in_edge_m[x])
								{
									break;
								}

							}
							if (x == in_edge_n.size())
							{
								in = true;
								if (out_edge_n.size() == out_edge_m.size())
								{
									int y = 0;
									for (y; y < out_edge_n.size(); y++)
									{
										if (out_edge_n[y] != out_edge_m[y])
										{
											break;
										}
									}
									if (y == out_edge_m.size()) out = true;
								}
							}
						}
						if (in && out) 
						{
							if (tran_seq_exist.count(pair1))
							{
								Cmp_Node& last_node = cmp_nodes.back();
								last_node.timestamps.push_back(tran_nodes[pair2].timestamp);
							}
							else
							{
								tran_seq_exist.insert(pair1);
								Cmp_Node new_node = Cmp_Node(cmp_nodes.size(), { tran_nodes[pair1].timestamp,tran_nodes[pair2].timestamp }, tran_nodes[pair1].identity, {});
								cmp_nodes.push_back(new_node);
								tran_seq_to_cmp_seq.insert(make_pair(tran_nodes[pair1].seq, new_node.seq));
							}
							next_n = m + 1;
						}
						else if (in)
						{
							if (tran_seq_exist.count(pair1))
							{
								tran_seq_exist.insert(pair2);
								Cmp_Node node1 = cmp_nodes.back();
								(*cmp_nodes.rbegin()).adj_edge.push_back(Cmp_Edge(node1.seq, node1.seq + 1));
								Cmp_Node node2 = Cmp_Node(cmp_nodes.size(), { tran_nodes[pair2].timestamp }, tran_nodes[pair2].identity, {});
								cmp_nodes.push_back(node2);
								tran_seq_to_cmp_seq.insert(make_pair(tran_nodes[pair2].seq, node2.seq));
								vertical_tran_seq.insert(tran_nodes[pair2].seq);
								int join_pos = 0;
								for (int k = out_edge_n.size() - 1; k >= 0; k--)
								{
									int p = out_edge_n[k];
									int l = out_edge_m.size() - 1;
									for (l; l >= 0; l--)
									{
										if (out_edge_m[l] == out_edge_n[k])
											break;
									}
									
									if (l == -1)
									{
										join_pos = k + 1;
										break;
									}
								}
								vector<int> tmp_out_edge;
								for (int k = 0; k < join_pos; k++)
								{
									tmp_out_edge.push_back(out_edge_n[k]);
								}
								out_edge_n = tmp_out_edge;
							}
							else
							{
								tran_seq_exist.insert(pair1);
								tran_seq_exist.insert(pair2);
								int size = cmp_nodes.size();
								Cmp_Node node1 = Cmp_Node(cmp_nodes.size(), { tran_nodes[pair1].timestamp }, tran_nodes[pair1].identity, {});
								node1.adj_edge.push_back(Cmp_Edge(node1.seq, node1.seq + 1));
								cmp_nodes.push_back(node1);
								Cmp_Node node2 = Cmp_Node(cmp_nodes.size(), { tran_nodes[pair2].timestamp }, tran_nodes[pair2].identity, {});
								cmp_nodes.push_back(node2);
								tran_seq_to_cmp_seq.insert(make_pair(tran_nodes[pair1].seq, node1.seq));
								tran_seq_to_cmp_seq.insert(make_pair(tran_nodes[pair2].seq, node2.seq));
								vertical_tran_seq.insert(tran_nodes[pair2].seq);
								int join_pos = 0;
								for (int k = out_edge_n.size() - 1; k >= 0; k--)
								{
									int p = out_edge_n[k];
									int l = out_edge_m.size() - 1;
									for (l; l >= 0; l--)
									{
										if (out_edge_m[l] == out_edge_n[k])
											break;
									}
									
									if (l == -1)
									{
										join_pos = k + 1;
										break;
									}
								}
								vector<int> tmp_out_edge;
								for (int k = 0; k < join_pos; k++)
								{
									tmp_out_edge.push_back(out_edge_n[k]);
								}
								out_edge_n = tmp_out_edge;
							}
							next_n = m;
							break;
						}
						else
						{
							if (tran_seq_exist.count(pair1))
							{
								tran_seq_exist.insert(pair2);
								Cmp_Node node2 = Cmp_Node(cmp_nodes.size(), { tran_nodes[pair2].timestamp }, tran_nodes[pair2].identity, {});
								cmp_nodes.push_back(node2);
								tran_seq_to_cmp_seq.insert(make_pair(tran_nodes[pair2].seq, node2.seq));
							}
							else
							{
								tran_seq_exist.insert(pair1);
								tran_seq_exist.insert(pair2);
								Cmp_Node node1 = Cmp_Node(cmp_nodes.size(), { tran_nodes[pair1].timestamp }, tran_nodes[pair1].identity, {});
								cmp_nodes.push_back(node1);
								Cmp_Node node2 = Cmp_Node(cmp_nodes.size(), { tran_nodes[pair2].timestamp }, tran_nodes[pair2].identity, {});
								cmp_nodes.push_back(node2);
								tran_seq_to_cmp_seq.insert(make_pair(tran_nodes[pair1].seq, node1.seq));
								tran_seq_to_cmp_seq.insert(make_pair(tran_nodes[pair2].seq, node2.seq));

							}
							next_n = m;
							break;
						}
					}
					if (n == next_n)
						break;
					else
						n = next_n;
				}
			}
		}
		for (auto& pair : tran_seq_to_cmp_seq)
		{
			int from_tran_seq = pair.first;
			int from_cmp_seq = pair.second;
			for (int j = 0; j < out_edges[tran_nodes[from_tran_seq].identity][from_tran_seq].size(); j++)
			{
				int to_tran_seq = out_edges[tran_nodes[from_tran_seq].identity][from_tran_seq][j];
				if (tran_seq_to_cmp_seq.count(to_tran_seq))
				{
					int to_cmp_seq = tran_seq_to_cmp_seq[to_tran_seq];
					
					if (vertical_tran_seq.count(to_tran_seq))
					{
						
					}
					else
					{
						Cmp_Node& cmp_node = cmp_nodes[from_cmp_seq];
						cmp_node.adj_edge.push_back(Cmp_Edge(from_cmp_seq, to_cmp_seq));
						cmp_edge_size++;
					}
				}
				else
				{

				}
			}
		}
		
		tran_graph.Tran_Data_Nodes.clear();

	}
};