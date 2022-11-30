#pragma once
#include"Transformer.h"
#include <stack>
#include <deque> 
#include <set>

using namespace std;
// set thread num
#define THREAD_NUM 1
struct  Computer
{
	Transformer transformer;
	vector<double> BC;
	TransformGraph& Tran_Graph = transformer.tran_graph;
	CompressGraph& Cmp_Graph = transformer.cmp_graph;
	unordered_map<string, int>& id_to_seq = transformer.data_graph.identity_to_seq;
	int data_nodes_size = transformer.data_graph.graph_data.size();
	unordered_map<int, int>& tran_root_seqs = transformer.tran_root_seqs;
	Computer(Transformer& t) :
		transformer(t) {}
	//compute shortest path and temporal betweenness with time instance graph
	void computeByShortestPathT()
	{
		 vector<Tran_Node>& tran_nodes = Tran_Graph.Tran_Data_Nodes;
		 int max_ts = transformer.data_graph.max_ts;
		 int tran_nodes_size = tran_nodes.size();
		 BC.clear();
		 BC.resize(data_nodes_size);
		 vector<vector<int>> t_distance(THREAD_NUM, vector<int>(data_nodes_size, -1));
		 vector<vector<int>> t_tran_distance(THREAD_NUM, vector<int>(tran_nodes_size, -1));
		 vector<vector<double>> t_tmp_sigma(THREAD_NUM, vector<double>(tran_nodes_size, 0));
		 vector<vector<int>> t_total_sigma(THREAD_NUM, vector<int>(tran_nodes_size, 0));
		 vector<vector<bool>> t_visited_tran_node(THREAD_NUM, vector<bool>(tran_nodes_size, 0));
		 vector<vector<vector<int>>> t_pre(THREAD_NUM, vector<vector<int>>(tran_nodes_size, vector<int>()));
		 omp_set_num_threads(THREAD_NUM);
		 #pragma omp parallel for
		 for (int i = 0; i < data_nodes_size; i++)
		 {
		 	int cur_thread_num = omp_get_thread_num();
		 	vector<int>& distance = t_distance[cur_thread_num];
		 	vector<int>& tran_distance = t_tran_distance[cur_thread_num];
		 	vector<double>& tmp_sigma = t_tmp_sigma[cur_thread_num];
		 	vector<vector<int>>& pre = t_pre[cur_thread_num];
		 	unordered_set<int> visited_tran_node ;
		 	transform(distance.cbegin(), distance.cend(), distance.begin(), [](auto l) { return -1; });
		 	transform(tran_distance.cbegin(), tran_distance.cend(), tran_distance.begin(), [](auto l) { return -1; });
		 	transform(tmp_sigma.cbegin(), tmp_sigma.cend(), tmp_sigma.begin(), [](auto l) { return 0; });
		 	transform(pre.cbegin(), pre.cend(), pre.begin(), [](auto l) { return vector<int>(); });
		 	//transform(visited_tran_node.cbegin(), visited_tran_node.cend(), visited_tran_node.begin(), [](auto l) { return 0; });
		 	deque<int> q;
		 	stack<int> stk;
		 	map<int, int> total_sigma;
		 	unordered_map<int, int> tran_sigma;
		 	unordered_map<int, double> tran_delta;
		 	Tran_Node& tran_root_node = tran_nodes[tran_root_seqs[i]];
		 	q.push_back(tran_root_node.seq);
		 	distance[tran_root_node.identity] = 0;
		 	tran_distance[tran_root_node.seq] = 0;
		 	tran_sigma[tran_root_node.seq] = 1;
		 	tmp_sigma[tran_root_node.seq] = 1;
		 	while (!q.empty())
		 	{
		 		int u = q.front();
		 		stk.push(u);
		 		q.pop_front();
		 		for (int j = 0; j < tran_nodes[u].adj_edge.size(); j++)
		 		{
		 			int u_ = tran_nodes[u].adj_edge[j].destination_seq;
		 			int t_ = tran_nodes[u_].timestamp;
		 			if (tran_distance[u_] == -1 || tran_distance[u] + 1 == tran_distance[u_])
		 			{
		 				tran_distance[u_] = tran_distance[u] + 1;
		 				tran_sigma[u_] += tmp_sigma[u];
		 				tmp_sigma[u_] += tmp_sigma[u];
		 				if (distance[tran_nodes[u_].identity] == -1)
		 				{
		 					distance[tran_nodes[u_].identity] = tran_distance[u_];
		 				}
		 				if (!visited_tran_node.count(u_))
		 				{
		 					q.push_back(u_);
		 					visited_tran_node.insert(u_);
		 				}
		 				pre[u_].push_back(u);
		 			}
		 		}
		 	}
		 	for (auto& s : tran_sigma)
		 	{
		 		if (tran_distance[s.first] == distance[tran_nodes[s.first].identity])
		 		{

		 			total_sigma[tran_nodes[s.first].identity] += s.second;
		 		}
		 	}
		 	while (!stk.empty())
		 	{
		 		int u_ = stk.top();
		 		stk.pop();
		 		for (int j = 0; j < pre[u_].size(); j++)
		 		{
		 			if (pre[u_][j] == tran_root_node.seq) continue;
		 			pair<int, int> ut = { pre[u_][j] , tran_nodes[pre[u_][j]].timestamp };
		 			if (tran_nodes[ut.first].identity == tran_nodes[u_].identity)
		 			{
		 				tran_delta[ut.first] += tran_delta[u_];
		 			}
		 			else
		 			{
		 				if (tran_distance[u_] == distance[tran_nodes[u_].identity])
		 				{
		 					tran_delta[ut.first] += ((double)(tran_sigma[ut.first]) / total_sigma[tran_nodes[u_].identity]);
		 				}
		 				tran_delta[ut.first] += (double)tran_delta[u_] * tran_sigma[ut.first] / tran_sigma[u_];
		 			}
		 		}
		 	}
		 	#pragma omp critical
		 	{
		 		for (auto& d : tran_delta)
		 		{
		 			BC[tran_nodes[d.first].identity] += d.second;
		 		}
		 	}
		 }
	}
	//compute shortest path and temporal betweenness by compression
	void computeByShortestPath()
	{
		vector<Cmp_Node>& cmp_nodes = Cmp_Graph.cmp_data_nodes;
		int max_ts = transformer.data_graph.max_ts;
		int cmp_nodes_size = cmp_nodes.size();
		BC.clear();
		BC.resize(data_nodes_size);
		vector<vector<int>> t_distance(THREAD_NUM, vector<int>(data_nodes_size, -1));
		vector<vector<int>> t_cmp_distance(THREAD_NUM, vector<int>(cmp_nodes_size, -1));
		vector<vector<int>> t_tmp_sigma(THREAD_NUM, vector<int>(cmp_nodes_size, 0));
		vector<vector<int>> t_total_sigma(THREAD_NUM, vector<int>(cmp_nodes.size(), 0));
		vector<vector<vector<int>>> t_pre(THREAD_NUM, vector<vector<int>>(cmp_nodes_size, vector<int>()));
		vector<vector<int>> t_visited_cmp_node(THREAD_NUM, vector<int>(cmp_nodes_size, 0));
		omp_set_num_threads(THREAD_NUM);
		#pragma omp parallel for
		for (int i = 0; i < data_nodes_size; i++)
		{
			int cur_thread_num = omp_get_thread_num();
			vector<int>& distance = t_distance[cur_thread_num];
			vector<int>& cmp_distance = t_cmp_distance[cur_thread_num];
			vector<int>& tmp_sigma = t_tmp_sigma[cur_thread_num];
			vector<vector<int>>& pre = t_pre[cur_thread_num];
			vector<int>& visited_cmp_node = t_visited_cmp_node[cur_thread_num];
			Cmp_Node& cmp_root_node = cmp_nodes[i];
			transform(distance.cbegin(), distance.cend(), distance.begin(), [](auto l) { return -1; });
			transform(cmp_distance.cbegin(), cmp_distance.cend(), cmp_distance.begin(), [](auto l) { return -1; });
			transform(tmp_sigma.cbegin(), tmp_sigma.cend(), tmp_sigma.begin(), [](auto l) { return 0; });
			transform(pre.cbegin(), pre.cend(), pre.begin(), [](auto l) { return vector<int>(); });
			transform(visited_cmp_node.cbegin(), visited_cmp_node.cend(), visited_cmp_node.begin(), [](auto l) { return 0; });
			unordered_map<int, int> total_sigma;
			unordered_map<int, int> cmp_sigma;
			unordered_map<int, double> cmp_delta;
			deque<int> q;
			stack<int> stk;
			q.push_back(cmp_root_node.seq);
			distance[cmp_root_node.identity] = 0;
			cmp_distance[cmp_root_node.seq] = 0;
			cmp_sigma[cmp_root_node.seq] = 1;
			tmp_sigma[cmp_root_node.seq] = 1;
			while (!q.empty())
			{
				int u = q.front();
				stk.push(u);
				q.pop_front();
				for (int j = 0; j < cmp_nodes[u].adj_edge.size(); j++)
				{
					int u_ = cmp_nodes[u].adj_edge[j].destination_seq;
					if (cmp_distance[u_] == -1 || cmp_distance[u] + 1 == cmp_distance[u_])
					{
						if (cmp_nodes[u].identity == cmp_nodes[u_].identity)
						{
							cmp_distance[u_] = cmp_distance[u];
							cmp_sigma[u_] += cmp_sigma[u];
							tmp_sigma[u_] += (cmp_sigma[u_] * cmp_nodes[u_].timestamps.size() + tmp_sigma[u]);
							if (!visited_cmp_node[u_])
							{
								q.push_front(u_);
								visited_cmp_node[u_] = 1;
							}
							for (int k = 0; k < pre[u].size(); k++)
							{
								if (cmp_nodes[pre[u][k]].identity != cmp_nodes[u].identity)
								{
									pre[u_].push_back(pre[u][k]);
								}
							}
						}
						else
						{
							cmp_distance[u_] = cmp_distance[u] + 1;
							cmp_sigma[u_] += tmp_sigma[u];
							tmp_sigma[u_] += tmp_sigma[u] * cmp_nodes[u_].timestamps.size();
							if (distance[cmp_nodes[u_].identity] == -1)
							{
								distance[cmp_nodes[u_].identity] = cmp_distance[u_];
							}
							if (!visited_cmp_node[u_])
							{
								q.push_back(u_);
								visited_cmp_node[u_] = 1;
							}
						}
						pre[u_].push_back(u);
					}
				}
			}
			for (auto& s : cmp_sigma)
			{
				if (cmp_distance[s.first] == distance[cmp_nodes[s.first].identity])
				{

					total_sigma[cmp_nodes[s.first].identity] += s.second * cmp_nodes[s.first].timestamps.size();
				}
			}
			while (!stk.empty())
			{
				int u_ = stk.top();
				stk.pop();
				for (int j = 0; j < pre[u_].size(); j++)
				{
					if (pre[u_][j] == cmp_root_node.seq) continue;
					pair<int, int> ut = { pre[u_][j] , cmp_nodes[pre[u_][j]].timestamps.front() };
					if (cmp_nodes[ut.first].identity == cmp_nodes[u_].identity)
					{
						cmp_delta[ut.first] += (cmp_delta[u_] * cmp_nodes[ut.first].timestamps.size() / cmp_nodes[u_].timestamps.size());
					}
					else
					{
						if (cmp_distance[u_] == distance[cmp_nodes[u_].identity])
						{
							cmp_delta[ut.first] += ((double)(cmp_sigma[ut.first] * cmp_nodes[u_].timestamps.size() * cmp_nodes[ut.first].timestamps.size()) / total_sigma[cmp_nodes[u_].identity]);
						}
						cmp_delta[ut.first] += (double)cmp_delta[u_] * cmp_sigma[ut.first] * cmp_nodes[ut.first].timestamps.size() / cmp_sigma[u_];
					}
				}
			}
			#pragma omp critical
			{
				for (auto& d : cmp_delta)
				{
					BC[cmp_nodes[d.first].identity] += d.second;
				}
			}
		}
	}
	//compute earliest path and temporal betweenness by compression
	void computeByForemostPath()
	{
		vector<Cmp_Node>& cmp_nodes = Cmp_Graph.cmp_data_nodes;
		int max_ts = transformer.data_graph.max_ts;
		int cmp_nodes_size = cmp_nodes.size();
		BC.clear();
		BC.resize(data_nodes_size);
		vector<vector<int>> t_fm_time(THREAD_NUM, vector<int>(data_nodes_size, 0));
		vector<vector<int>> t_tmp_sigma(THREAD_NUM, vector<int>(cmp_nodes_size, 0));
		vector<vector<int>> t_total_sigma(THREAD_NUM, vector<int>(cmp_nodes.size(), 0));
		vector<vector<int>> t_bfs(THREAD_NUM, vector<int>(cmp_nodes.size(), -1));
		vector<vector<vector<int>>> t_pre(THREAD_NUM, vector<vector<int>>(cmp_nodes_size, vector<int>()));
		omp_set_num_threads(THREAD_NUM);
		#pragma omp parallel for
		for (int i = 0; i < data_nodes_size; i++)
		{
			int cur_thread_num = omp_get_thread_num();
			Cmp_Node& cmp_root_node = cmp_nodes[i];
			vector<int>& fm_time = t_fm_time[cur_thread_num];
			transform(fm_time.cbegin(), fm_time.cend(), fm_time.begin(), [](auto l) { return INT_MAX; });
			unordered_set<int> visited_cmp;
			queue<int> q;
			for (int j = 0; j < cmp_root_node.adj_edge.size(); j++)
			{
				int des = cmp_root_node.adj_edge[j].destination_seq;
				q.push(des);
				visited_cmp.insert(des);
				fm_time[cmp_nodes[des].identity] = min(fm_time[cmp_nodes[des].identity], cmp_nodes[des].timestamps[0]);
			}
			while (!q.empty())
			{
				int u = q.front();
				q.pop();
				for (int j = 0; j < cmp_nodes[u].adj_edge.size(); j++)
				{
					int des = cmp_nodes[u].adj_edge[j].destination_seq;
					if (!visited_cmp.count(des) && cmp_nodes[u].identity != cmp_nodes[des].identity)
					{
						q.push(des);
						visited_cmp.insert(des);
						fm_time[cmp_nodes[des].identity] = min(fm_time[cmp_nodes[des].identity], cmp_nodes[des].timestamps[0]);
					}
				}
			}
			deque<pair<int, int>> dq;
			stack<int> stk;
			vector<int>& total_sigma = t_total_sigma[cur_thread_num];
			unordered_map<int, int> cmp_sigma;
			vector<int>& tmp_sigma = t_tmp_sigma[cur_thread_num];
			unordered_map<int, int> update_tmp_sigma;
			unordered_map<int, int> update_cmp_sigma;
			vector<vector<int>>& pre = t_pre[cur_thread_num];
			vector<int>& bfs = t_bfs[cur_thread_num];
			unordered_map<int, double> cmp_delta;
			transform(total_sigma.cbegin(), total_sigma.cend(), total_sigma.begin(), [](auto l) { return 0; });
			transform(tmp_sigma.cbegin(), tmp_sigma.cend(), tmp_sigma.begin(), [](auto l) { return 0; });
			transform(bfs.cbegin(), bfs.cend(), bfs.begin(), [](auto l) { return -1; });
			transform(pre.cbegin(), pre.cend(), pre.begin(), [](auto l) { return vector<int>(); });
			visited_cmp.clear();
			dq.push_back({ cmp_root_node.seq ,0 });
			cmp_sigma.insert({ cmp_root_node.seq,1 });
			tmp_sigma[cmp_root_node.seq] = 1;
			bfs[cmp_root_node.seq] = 0;
			int last = -1;
			while (!dq.empty())
			{
				pair<int, int> u = dq.front();
				if (last == u.first)
				{
					dq.pop_front();
					continue;
				}
				last = u.first;
				stk.push(u.first);
				dq.pop_front();
				for (int j = 0; j < cmp_nodes[u.first].adj_edge.size(); j++)
				{
					int u_ = cmp_nodes[u.first].adj_edge[j].destination_seq;
					if (u.second == 0)
					{
						if (bfs[u_] != -1 && bfs[u_] <= bfs[u.first])
						{
							cmp_sigma[u_] += tmp_sigma[u.first];
							update_cmp_sigma[u_] += tmp_sigma[u.first];
							tmp_sigma[u_] += tmp_sigma[u.first] * cmp_nodes[u_].timestamps.size();
							update_tmp_sigma[u_] += tmp_sigma[u.first] * cmp_nodes[u_].timestamps.size();
							dq.push_back({ u_ ,1 });
							if (fm_time[cmp_nodes[u_].identity] == cmp_nodes[u_].timestamps[0])
							{
								total_sigma[cmp_nodes[u_].identity] += tmp_sigma[u.first];
							}
							pre[u_].push_back(u.first);
						}
						else
						{
							if (cmp_nodes[u.first].identity == cmp_nodes[u_].identity && (bfs[u_] == -1 || bfs[u.first] == bfs[u_]))
							{
								bfs[u_] = bfs[u.first];
								cmp_sigma[u_] += cmp_sigma[u.first];
								tmp_sigma[u_] += cmp_sigma[u_] * cmp_nodes[u_].timestamps.size() + tmp_sigma[u.first];
								if (!visited_cmp.count(u_))
								{
									visited_cmp.insert(u_);
									dq.push_front({ u_ ,0 });
								}
								for (int k = 0; k < pre[u.first].size(); k++)
								{
									pre[u_].push_back(pre[u.first][k]);
								}
							}
							else if (cmp_nodes[u.first].identity != cmp_nodes[u_].identity && (bfs[u_] == -1 || bfs[u.first] + 1 == bfs[u_]))
							{
								bfs[u_] = bfs[u.first] + 1;
								cmp_sigma[u_] += tmp_sigma[u.first];
								tmp_sigma[u_] += tmp_sigma[u.first] * cmp_nodes[u_].timestamps.size();
								if (!visited_cmp.count(u_))
								{
									visited_cmp.insert(u_);
									dq.push_back({ u_ ,0 });
								}
								if (fm_time[cmp_nodes[u_].identity] == cmp_nodes[u_].timestamps[0])
								{
									total_sigma[cmp_nodes[u_].identity] += tmp_sigma[u.first];
								}
							}
							pre[u_].push_back(u.first);
						}
					}
					else
					{
						if (cmp_nodes[u.first].identity == cmp_nodes[u_].identity && (bfs[u_] == -1 || bfs[u.first] == bfs[u_]))
						{
							cmp_sigma[u_] += update_cmp_sigma[u.first];
							update_cmp_sigma[u_] = update_cmp_sigma[u.first];
							tmp_sigma[u_] += update_cmp_sigma[u_] * cmp_nodes[u_].timestamps.size() + update_tmp_sigma[u.first];
							update_tmp_sigma[u_] = update_cmp_sigma[u_] * cmp_nodes[u_].timestamps.size() + update_tmp_sigma[u.first];
							dq.push_front({ u_ ,1 });
						}
						else if (cmp_nodes[u.first].identity != cmp_nodes[u_].identity && (bfs[u_] == -1 || bfs[u.first] + 1 == bfs[u_]))
						{
							cmp_sigma[u_] += update_tmp_sigma[u.first];
							update_cmp_sigma[u_] += update_tmp_sigma[u.first];
							tmp_sigma[u_] += update_tmp_sigma[u.first] * cmp_nodes[u_].timestamps.size();
							update_tmp_sigma[u_] += update_tmp_sigma[u.first] * cmp_nodes[u_].timestamps.size();
							if (fm_time[cmp_nodes[u_].identity] == cmp_nodes[u_].timestamps[0])
							{
								total_sigma[cmp_nodes[u_].identity] += update_tmp_sigma[u.first];
							}
							dq.push_back({ u_ ,1 });
						}
					}
					if (j == cmp_nodes[u.first].adj_edge.size() - 1)
					{
						update_cmp_sigma[u.first] = 0;
						update_tmp_sigma[u.first] = 0;
					}
				}
			}
			visited_cmp.clear();
			while (!stk.empty())
			{
				int u_ = stk.top();
				stk.pop();
				if (visited_cmp.count(u_)) continue;
				visited_cmp.insert(u_);
				bool flag = fm_time[cmp_nodes[u_].identity] == cmp_nodes[u_].timestamps[0];
				for (int j = 0; j < pre[u_].size(); j++)
				{
					if (pre[u_][j] == cmp_root_node.seq) continue;
					int u = pre[u_][j];
					if (cmp_nodes[u].identity == cmp_nodes[u_].identity)
					{
						cmp_delta[u] += (cmp_delta[u_] * cmp_nodes[u].timestamps.size() / cmp_nodes[u_].timestamps.size());
					}
					else
					{
						if (flag)
						{
							cmp_delta[u] += ((double)(cmp_sigma[u] * cmp_nodes[u].timestamps.size()) / total_sigma[cmp_nodes[u_].identity]);
						}
						cmp_delta[u] += (double)cmp_delta[u_] * cmp_sigma[u] * cmp_nodes[u].timestamps.size() / cmp_sigma[u_];
					}
				}
			}
#pragma omp critical
			{
				for (auto& d : cmp_delta)
				{
					BC[cmp_nodes[d.first].identity] += d.second;
				}
			}
		}
	}
	//compute earliest path and temporal betweenness with time instance graph
	void computeByForemostPathT()
	{
		vector<Tran_Node>& tran_nodes = Tran_Graph.Tran_Data_Nodes;
		int max_ts = transformer.data_graph.max_ts;
		int tran_nodes_size = tran_nodes.size();
		BC.clear();
		BC.resize(data_nodes_size);
		vector<vector<int>> t_fm_time(THREAD_NUM, vector<int>(data_nodes_size, 0));
		vector<vector<int>> t_tmp_sigma(THREAD_NUM, vector<int>(tran_nodes_size, 0));
		vector<vector<int>> t_total_sigma(THREAD_NUM, vector<int>(tran_nodes_size, 0));
		vector<vector<int>> t_bfs(THREAD_NUM, vector<int>(tran_nodes_size, -1));
		vector<vector<vector<int>>> t_pre(THREAD_NUM, vector<vector<int>>(tran_nodes_size, vector<int>()));
		omp_set_num_threads(THREAD_NUM);
#pragma omp parallel for
		for (int i = 0; i < data_nodes_size; i++)
		{
			int cur_thread_num = omp_get_thread_num();
			Tran_Node& tran_root_node = tran_nodes[tran_root_seqs[i]];
			vector<int>& fm_time = t_fm_time[cur_thread_num];
			transform(fm_time.cbegin(), fm_time.cend(), fm_time.begin(), [](auto l) { return INT_MAX; });
			unordered_set<int> visited_tran;
			queue<int> q;
			for (int j = 0; j < tran_root_node.adj_edge.size(); j++)
			{
				int des = tran_root_node.adj_edge[j].destination_seq;
				q.push(des);
				visited_tran.insert(des);
				fm_time[tran_nodes[des].identity] = min(fm_time[tran_nodes[des].identity], tran_nodes[des].timestamp);
			}
			while (!q.empty())
			{
				int u = q.front();
				q.pop();
				for (int j = 0; j < tran_nodes[u].adj_edge.size(); j++)
				{
					int des = tran_nodes[u].adj_edge[j].destination_seq;
					if (!visited_tran.count(des) && tran_nodes[u].identity != tran_nodes[des].identity)
					{
						q.push(des);
						visited_tran.insert(des);
						fm_time[tran_nodes[des].identity] = min(fm_time[tran_nodes[des].identity], tran_nodes[des].timestamp);
					}
				}
			}
			deque<pair<int, int>> dq;//
			stack<int> stk;
			vector<int>& total_sigma = t_total_sigma[cur_thread_num];
			unordered_map<int, int> tran_sigma;
			vector<int>& tmp_sigma = t_tmp_sigma[cur_thread_num];
			unordered_map<int, int> update_tmp_sigma;
			unordered_map<int, int> update_tran_sigma;
			vector<vector<int>>& pre = t_pre[cur_thread_num];
			vector<int>& bfs = t_bfs[cur_thread_num];
			unordered_map<int, double> tran_delta;
			transform(total_sigma.cbegin(), total_sigma.cend(), total_sigma.begin(), [](auto l) { return 0; });
			transform(tmp_sigma.cbegin(), tmp_sigma.cend(), tmp_sigma.begin(), [](auto l) { return 0; });
			transform(bfs.cbegin(), bfs.cend(), bfs.begin(), [](auto l) { return -1; });
			transform(pre.cbegin(), pre.cend(), pre.begin(), [](auto l) { return vector<int>(); });
			visited_tran.clear();
			dq.push_back({ tran_root_node.seq ,0 });
			tran_sigma.insert({ tran_root_node.seq,1 });
			tmp_sigma[tran_root_node.seq] = 1;
			bfs[tran_root_node.seq] = 0;
			int last = -1;
			while (!dq.empty())
			{
				pair<int, int> u = dq.front();
				if (last == u.first)
				{
					dq.pop_front();
					continue;
				}
				last = u.first;
				stk.push(u.first);
				dq.pop_front();
				for (int j = 0; j < tran_nodes[u.first].adj_edge.size(); j++)
				{
					int u_ = tran_nodes[u.first].adj_edge[j].destination_seq;
					if (u.second == 0)
					{
						if (bfs[u_] != -1 && bfs[u_] <= bfs[u.first])
						{
							tran_sigma[u_] += tmp_sigma[u.first];
							update_tran_sigma[u_] += tmp_sigma[u.first];
							tmp_sigma[u_] += tmp_sigma[u.first];
							update_tmp_sigma[u_] += tmp_sigma[u.first];
							dq.push_back({ u_ ,1 });
							if (fm_time[tran_nodes[u_].identity] == tran_nodes[u_].timestamp)
							{
								total_sigma[tran_nodes[u_].identity] += tmp_sigma[u.first];
							}
							pre[u_].push_back(u.first);
						}
						else
						{
							bfs[u_] = bfs[u.first] + 1;
							tran_sigma[u_] += tmp_sigma[u.first];
							tmp_sigma[u_] += tmp_sigma[u.first];
							if (!visited_tran.count(u_))
							{
								visited_tran.insert(u_);
								dq.push_back({ u_ ,0 });
							}
							if (fm_time[tran_nodes[u_].identity] == tran_nodes[u_].timestamp)
							{
								total_sigma[tran_nodes[u_].identity] += tmp_sigma[u.first];
							}
							pre[u_].push_back(u.first);
						}
					}
					else
					{
						if (bfs[u_] == -1 || bfs[u.first] + 1 == bfs[u_])
						{
							tran_sigma[u_] += update_tmp_sigma[u.first];
							update_tran_sigma[u_] += update_tmp_sigma[u.first];
							tmp_sigma[u_] += update_tmp_sigma[u.first];
							update_tmp_sigma[u_] += update_tmp_sigma[u.first];
							if (fm_time[tran_nodes[u_].identity] == tran_nodes[u_].timestamp)
							{
								total_sigma[tran_nodes[u_].identity] += update_tmp_sigma[u.first];
							}
							dq.push_back({ u_ ,1 });
						}

					}
				}
				update_tran_sigma[u.first] = 0;
				update_tmp_sigma[u.first] = 0;
			}
			visited_tran.clear();
			while (!stk.empty())
			{
				int u_ = stk.top();
				stk.pop();
				if (visited_tran.count(u_)) continue;
				visited_tran.insert(u_);
				bool flag = fm_time[tran_nodes[u_].identity] == tran_nodes[u_].timestamp;
				for (int j = 0; j < pre[u_].size(); j++)
				{
					if (pre[u_][j] == tran_root_node.seq) continue;
					int u = pre[u_][j];
					if (flag)
					{
						tran_delta[u] += ((double)(tran_sigma[u]) / total_sigma[tran_nodes[u_].identity]);
					}
					tran_delta[u] += (double)tran_delta[u_] * tran_sigma[u] / tran_sigma[u_];
				}
			}
#pragma omp critical
			{
				for (auto& d : tran_delta)
				{
					BC[tran_nodes[d.first].identity] += d.second;
				}
			}
		}
	}
	//compute shortest and earliest path and temporal betweenness by compression
	void computeByShortestForemostPath()
	{
		vector<Cmp_Node>& cmp_nodes = Cmp_Graph.cmp_data_nodes;
		int max_ts = transformer.data_graph.max_ts;
		int cmp_nodes_size = cmp_nodes.size();
		BC.resize(data_nodes_size);
		vector<vector<int>> t_distance(THREAD_NUM, vector<int>(data_nodes_size, -1));
		vector<vector<int>> t_cmp_distance(THREAD_NUM, vector<int>(cmp_nodes_size, -1));
		vector<vector<double>> t_tmp_sigma(THREAD_NUM, vector<double>(cmp_nodes_size, 0));
		vector<vector<double>> t_total_delta(THREAD_NUM, vector<double>(cmp_nodes.size(), 0));
		vector<vector<int>> t_total_sigma(THREAD_NUM, vector<int>(cmp_nodes.size(), 0));
		vector<vector<vector<int>>> t_pre(THREAD_NUM, vector<vector<int>>(cmp_nodes_size, vector<int>()));
		vector<vector<int>> t_fm_time(THREAD_NUM, vector<int>(data_nodes_size, INT_MAX));
		omp_set_num_threads(THREAD_NUM);
#pragma omp parallel for
		for (int i = 0; i < data_nodes_size; i++)
		{
			int cur_thread_num = omp_get_thread_num();
			vector<int>& distance = t_distance[cur_thread_num];
			vector<int>& cmp_distance = t_cmp_distance[cur_thread_num];
			vector<double>& tmp_sigma = t_tmp_sigma[cur_thread_num];
			vector<double>& total_delta = t_total_delta[cur_thread_num];
			vector<vector<int>>& pre = t_pre[cur_thread_num];
			vector<int>& fm_time = t_fm_time[cur_thread_num];
			Cmp_Node& cmp_root_node = cmp_nodes[i];
			transform(distance.cbegin(), distance.cend(), distance.begin(), [](auto l) { return -1; });
			transform(cmp_distance.cbegin(), cmp_distance.cend(), cmp_distance.begin(), [](auto l) { return -1; });
			transform(tmp_sigma.cbegin(), tmp_sigma.cend(), tmp_sigma.begin(), [](auto l) { return 0; });
			transform(total_delta.cbegin(), total_delta.cend(), total_delta.begin(), [](auto l) { return 0; });
			transform(pre.cbegin(), pre.cend(), pre.begin(), [](auto l) { return vector<int>(); });
			transform(fm_time.cbegin(), fm_time.cend(), fm_time.begin(), [](auto l) { return INT_MAX; });
			unordered_set<int> visited_cmp_node;
			unordered_map<int, int> cmp_sigma;
			unordered_map<int, double> cmp_delta;
			deque<int> q;
			stack<int> stk;
			q.push_back(cmp_root_node.seq);
			distance[cmp_root_node.identity] = 0;
			cmp_distance[cmp_root_node.seq] = 0;
			cmp_sigma[cmp_root_node.seq] = 1;
			tmp_sigma[cmp_root_node.seq] = 1;
			while (!q.empty())
			{
				int u = q.front();
				stk.push(u);
				fm_time[cmp_nodes[u].identity] = min(fm_time[cmp_nodes[u].identity], cmp_nodes[u].timestamps[0]);
				q.pop_front();
				for (int j = 0; j < cmp_nodes[u].adj_edge.size(); j++)
				{
					int u_ = cmp_nodes[u].adj_edge[j].destination_seq;
					if (cmp_distance[u_] == -1 || cmp_distance[u] + 1 == cmp_distance[u_])
					{
						if (cmp_nodes[u].identity == cmp_nodes[u_].identity)
						{
							cmp_distance[u_] = cmp_distance[u];
							cmp_sigma[u_] += cmp_sigma[u];
							tmp_sigma[u_] += (cmp_sigma[u_] * cmp_nodes[u_].timestamps.size() + tmp_sigma[u]);
							if (!visited_cmp_node.count(u_))
							{
								q.push_front(u_);
								visited_cmp_node.insert(u_);
							}
							for (int k = 0; k < pre[u].size(); k++)
							{
								if (cmp_nodes[pre[u][k]].identity != cmp_nodes[u].identity)
								{
									pre[u_].push_back(pre[u][k]);
								}
							}
						}
						else
						{
							cmp_distance[u_] = cmp_distance[u] + 1;
							cmp_sigma[u_] += tmp_sigma[u];
							tmp_sigma[u_] += tmp_sigma[u] * cmp_nodes[u_].timestamps.size();
							if (distance[cmp_nodes[u_].identity] == -1)
							{
								distance[cmp_nodes[u_].identity] = cmp_distance[u_];
							}
							if (!visited_cmp_node.count(u_))
							{
								q.push_back(u_);
								visited_cmp_node.insert(u_);
							}
						}
						pre[u_].push_back(u);
					}
				}
			}
			while (!stk.empty())
			{
				int u_ = stk.top();
				stk.pop();
				for (int j = 0; j < pre[u_].size(); j++)
				{
					if (pre[u_][j] == cmp_root_node.seq) continue;
					pair<int, int> ut = { pre[u_][j] , cmp_nodes[pre[u_][j]].timestamps.front() };
					int flag = fm_time[cmp_nodes[u_].identity] == cmp_nodes[u_].timestamps[0];
					if (cmp_nodes[ut.first].identity == cmp_nodes[u_].identity)
					{
						cmp_delta[ut.first] += (cmp_delta[u_] * cmp_nodes[ut.first].timestamps.size() / cmp_nodes[u_].timestamps.size());
					}
					else
					{
						if (flag)
						{
							cmp_delta[ut.first] += ((double)(cmp_sigma[ut.first] * cmp_nodes[ut.first].timestamps.size()) / cmp_sigma[u_]);
						}
						cmp_delta[ut.first] += (double)cmp_delta[u_] * cmp_sigma[ut.first] * cmp_nodes[ut.first].timestamps.size() / cmp_sigma[u_];
					}
				}
			}
			#pragma omp critical
			{
				for (auto& d : cmp_delta)
				{
					BC[cmp_nodes[d.first].identity] += d.second;
				}
			}
		}
	}
	//compute shortest and earliest path and temporal betweenness with time instance graph
	void computeByShortestForemostPathT()
	{
		vector<Tran_Node>& tran_nodes = Tran_Graph.Tran_Data_Nodes;
		int max_ts = transformer.data_graph.max_ts;
		int tran_nodes_size = tran_nodes.size();
		BC.resize(data_nodes_size);
		vector<vector<int>> t_distance(THREAD_NUM, vector<int>(data_nodes_size, -1));
		vector<vector<int>> t_tran_distance(THREAD_NUM, vector<int>(tran_nodes_size, -1));
		vector<vector<double>> t_tmp_sigma(THREAD_NUM, vector<double>(tran_nodes_size, 0));
		vector<vector<vector<int>>> t_pre(THREAD_NUM, vector<vector<int>>(tran_nodes_size, vector<int>()));
		vector<vector<int>> t_fm_time(THREAD_NUM, vector<int>(data_nodes_size, INT_MAX));
		omp_set_num_threads(THREAD_NUM);
		#pragma omp parallel for
		for (int i = 0; i < data_nodes_size; i++)
		{
			int cur_thread_num = omp_get_thread_num();
			vector<int>& distance = t_distance[cur_thread_num];
			vector<int>& tran_distance = t_tran_distance[cur_thread_num];
			vector<double>& tmp_sigma = t_tmp_sigma[cur_thread_num];
			vector<vector<int>>& pre = t_pre[cur_thread_num];
			vector<int>& fm_time = t_fm_time[cur_thread_num];
			transform(distance.cbegin(), distance.cend(), distance.begin(), [](auto l) { return -1; });
			transform(tran_distance.cbegin(), tran_distance.cend(), tran_distance.begin(), [](auto l) { return -1; });
			transform(tmp_sigma.cbegin(), tmp_sigma.cend(), tmp_sigma.begin(), [](auto l) { return 0; });
			transform(pre.cbegin(), pre.cend(), pre.begin(), [](auto l) { return vector<int>(); });
			transform(fm_time.cbegin(), fm_time.cend(), fm_time.begin(), [](auto l) { return INT_MAX; });
			deque<int> q;
			stack<int> stk;
			unordered_map<int, int> total_sigma;
			unordered_set<int> visited_cmp_node;
			unordered_map<int, int> tran_sigma;
			unordered_map<int, double> tran_delta;
			Tran_Node& tran_root_node = tran_nodes[tran_root_seqs[i]];
			q.push_back(tran_root_node.seq);
			distance[tran_root_node.identity] = 0;
			tran_distance[tran_root_node.seq] = 0;
			tran_sigma[tran_root_node.seq] = 1;
			tmp_sigma[tran_root_node.seq] = 1;
			while (!q.empty())
			{
				int u = q.front();
				stk.push(u);
				fm_time[tran_nodes[u].identity] = min(fm_time[tran_nodes[u].identity], tran_nodes[u].timestamp);
				q.pop_front();
				for (int j = 0; j < tran_nodes[u].adj_edge.size(); j++)
				{
					int u_ = tran_nodes[u].adj_edge[j].destination_seq;
					int t_ = tran_nodes[u_].timestamp;
					if (tran_distance[u_] == -1 || tran_distance[u] + 1 == tran_distance[u_])
					{
						tran_distance[u_] = tran_distance[u] + 1;
						tran_sigma[u_] += tmp_sigma[u];
						tmp_sigma[u_] += tmp_sigma[u];
						if (distance[tran_nodes[u_].identity] == -1)
						{
							distance[tran_nodes[u_].identity] = tran_distance[u_];
						}

						if (!visited_cmp_node.count(u_))
						{
							q.push_back(u_);
							visited_cmp_node.insert(u_);
						}
						pre[u_].push_back(u);
					}
				}
			}
			while (!stk.empty())
			{
				int u_ = stk.top();
				stk.pop();
				for (int j = 0; j < pre[u_].size(); j++)
				{
					if (pre[u_][j] == tran_root_node.seq) continue;
					pair<int, int> ut = { pre[u_][j] , tran_nodes[pre[u_][j]].timestamp };
					int flag = fm_time[tran_nodes[u_].identity] == tran_nodes[u_].timestamp;
					if (flag)
					{
						tran_delta[ut.first] += ((double)(tran_sigma[ut.first]) / tran_sigma[u_]);
					}
					tran_delta[ut.first] += (double)tran_delta[u_] * tran_sigma[ut.first] / tran_sigma[u_];
				}
			}
			#pragma omp critical
			{
				for (auto& d : tran_delta)
				{
					BC[tran_nodes[d.first].identity] += d.second;
				}
			}
		}
	}

};