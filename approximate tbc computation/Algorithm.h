#pragma once

#include <iostream>
#include <utility>
#include <tuple>
#include <stack>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include <random>
#include <string>
#include <algorithm>

#include <numeric>
#include <cmath>
#include <ttmath/ttmath.hpp>

#include <omp.h>
#include <chrono>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/RadeAux.hpp>
#include <malloc.h>
#include "Transformer.h"

#define THREAD_NUM 24

using namespace std;
using namespace chrono;

struct Settings
{
	bool needCompress = true;
	bool isStrict = false;
    const double penaltyFactor = 1.0;
    const bool useRefinedEstimator=false;
    double epsilon = 0.01;
    double delta = 0.1;
};

struct Dsu {
    int n;
    vector<int> pa;
    Dsu(int n_ = 0) {
        n = n_;
        pa.resize(n);
        for (int i = 0; i < n; i++) {
            pa[i] = i;
        }
    }
    inline int find(int x) {
        return (pa[x] == x ? x : pa[x] = find(pa[x]));
    }
    inline bool unite(int x, int y) {
        int xr = find(x), yr = find(y);
        if (xr != yr) {
            pa[yr] = xr;
            return true;
        }
        return false;
    }
};

void get_topo(int data_nodes_size, Transformer &transformer, vector<pair<int, int>> &topo_seq) {
    /*
        get topological order
    */
    int cmp_edge_size = transformer.cmp_edge_size;
    vector<Cmp_Node> &cmp_nodes = transformer.cmp_graph.cmp_data_nodes;
    
    vector<int> data_in_edges(data_nodes_size, 0);
    for(auto &e : transformer.data_graph.edge_set) {
        data_in_edges[e.first.second]++;
    }
    queue<int> q;
    int sz = 0, cursor = data_nodes_size;
    vector<int> cmp_in_edges(cmp_nodes.size(), 0);
    for(int i=0; i<cmp_nodes.size(); ++i) {
        for(auto &e : cmp_nodes[i].adj_edge) {
            cmp_in_edges[e.destination_seq]++;
        }
    }
    for(int i=0; i<cmp_nodes.size(); ++i) {
        if(cmp_in_edges[i] == 0) {
            q.push(i);
        }
    }
    assert(!q.empty());
    int cnt = 0;
    while(!q.empty()) {
        auto cur = q.front();
        q.pop();
        int d_seq = cmp_nodes[cur].identity;
        if(topo_seq[d_seq].first == -1) {
            topo_seq[d_seq].first = cnt++;
        } else {
            topo_seq[d_seq].second = cnt++;
        }
        for(auto &e : cmp_nodes[cur].adj_edge) {
            cmp_in_edges[e.destination_seq]--;
            if(cmp_in_edges[e.destination_seq] == 0) {
                q.push(e.destination_seq);
            }
        }
    }
}


class RunningObjFuncSumExponents {
    private:
        map<size_t,double> objFuncSumExponentsMap;
    public:
        RunningObjFuncSumExponents() {}
        RunningObjFuncSumExponents(const map<size_t, double> &
                objFuncSumExponentsMap):
            objFuncSumExponentsMap(objFuncSumExponentsMap) {}

        RunningObjFuncSumExponents& operator+=(const RunningObjFuncSumExponents& rhs) {
            objFuncSumExponentsMap.insert(rhs.objFuncSumExponentsMap.begin(),
                    rhs.objFuncSumExponentsMap.end());
            return *this;
        }

        friend RunningObjFuncSumExponents operator+(RunningObjFuncSumExponents
                lhs, const RunningObjFuncSumExponents& rhs) {
            lhs += rhs;
            return lhs;
        }

        RunningObjFuncSumExponents& operator+=(const pair<size_t, double> &
                rhs) {
            objFuncSumExponentsMap.insert(rhs);
            return *this;
        }

        vector<double> getVector() {
            vector<double> v;
            for (auto entry : objFuncSumExponentsMap) {
                v.push_back(entry.second);
            }
            return v;
        }
};

uint64_t randomNode_temporal(uint64_t size) {
    auto &gen = Aux::Random::getURNG();
    uniform_int_distribution<uint64_t> distr{0, size - 1};  // random sample in range of [a, b]
    return distr(gen);
}

auto time_cost = [](const auto st, const auto en) {
	return ((double)duration_cast<microseconds>(en - st).count()) * microseconds::period::num / microseconds::period::den;
};

vector<int> get_topological_sequence(Transformer &transformer, Settings &setting)
{
    CompressGraph& Cmp_Graph = transformer.cmp_graph;
    vector<Cmp_Node>& cmp_nodes = Cmp_Graph.cmp_data_nodes;
    unordered_map<string, int>& id_to_seq = transformer.data_graph.identity_to_seq;
    unordered_map<int, string>& seq_to_id = transformer.data_graph.seq_to_identity;
    
    int data_nodes_size = transformer.data_graph.graph_data.size();		// size of t_nodes
    int cmp_nodes_size = cmp_nodes.size();
    vector<int> indegree(cmp_nodes_size, 0);
    for(int i=data_nodes_size; i<cmp_nodes_size; ++i)
    {
        for(auto& edge : cmp_nodes[i].adj_edge)
        {
            indegree[edge.destination_seq]++;
        }
    }
    vector<int> t_seq;
    queue<int> q;
    for(int i=data_nodes_size; i<cmp_nodes_size; ++i)
    {
        if(indegree[i] == 0)
        {
            q.push(i);
        }
    }
    int cur;
    while(!q.empty())
    {
        cur = q.front();
        q.pop();
        t_seq.push_back(cur);
        for(auto& edge : cmp_nodes[cur].adj_edge)
        {
            if((--indegree[edge.destination_seq]) == 0)
            {
                q.push(edge.destination_seq);
            }
        }
    }
    vector<int> v_seq;
    vector<bool> judged(data_nodes_size, false);
    for(auto& seq : t_seq)
    {
        cur = cmp_nodes[seq].identity;
        if(!judged[cur])
        {
            judged[cur] = true;
            v_seq.push_back(cur);
        }
    }
    return v_seq;
}



vector<double> sample_compress_v_v_parallel_purning_experiment(Transformer &transformer, Settings &setting, string dataName, Dsu &dsu)
{
    /*
        to run experiments and fix some data structure
    */
    omp_set_num_threads(THREAD_NUM);
    CompressGraph& Cmp_Graph = transformer.cmp_graph;
    vector<Cmp_Node>& cmp_nodes = Cmp_Graph.cmp_data_nodes;
    unordered_map<string, int>& id_to_seq = transformer.data_graph.identity_to_seq;
    unordered_map<int, string>& seq_to_id = transformer.data_graph.seq_to_identity;
    int max_ts = transformer.data_graph.max_ts;
    int cmp_nodes_size = cmp_nodes.size();
    int data_nodes_size = transformer.data_graph.graph_data.size();		// size of t_nodes
    int n = data_nodes_size;

    vector<double> scoreData(data_nodes_size, 0);
    vector<double> bcResult(data_nodes_size, 0);
    pair<double,double> minimizationResults(0, 100); // minimum and point of minimum for the objective function used in stopping condition.

    // data structure of bc calculation
    vector<double> sumContribsSquared(data_nodes_size, 0);
    vector<double> sumContribs(data_nodes_size, 0);
    vector<int> sumContribsCounts(data_nodes_size, 0);

    double supDeviationBound = 1;

    // sample number
    uint64_t sampleSize = ceil(log(2.0 / setting.delta) * (sqrt(1 + 16 * setting.epsilon)  + 1 + 8 * setting.epsilon) / (4 * setting.epsilon * setting.epsilon));
    uint64_t toSample = sampleSize;  // remember to fix!!!!!
    uint64_t iterations = 0;
    
    vector<vector<int>> t_distance(THREAD_NUM,vector<int>(data_nodes_size, -1));
    vector<vector<int>> t_cmp_distance(THREAD_NUM,vector<int>(cmp_nodes_size, -1));
    vector<vector<float>> t_tmp_sigma(THREAD_NUM,vector<float>(cmp_nodes_size, 0));
    vector<vector<float>> t_total_delta(THREAD_NUM,vector<float>(cmp_nodes.size(), 0));
    vector<vector<int>> t_total_sigma(THREAD_NUM, vector<int>(cmp_nodes.size(), 0));
    vector<vector<vector<int>>> t_pre(THREAD_NUM,vector<vector<int>>(cmp_nodes_size, vector<int>()));
    vector<unordered_map<int, int>> tt_total_sigma(THREAD_NUM);
    vector<unordered_set<int>> t_visited_cmp_node(THREAD_NUM);
    vector<unordered_map<int, int>> t_cmp_sigma(THREAD_NUM);
    vector<unordered_map<int, unordered_map<int, float>>> t_node_contribs(THREAD_NUM);
    double total_calculation_cost_time = 0;
    double total_rademacher_cost_time = 0;

    auto topo_start_time = chrono::steady_clock::now();
    vector<pair<int, int>> topo_seq(data_nodes_size, {-1, INT32_MAX});
    get_topo(data_nodes_size, transformer, topo_seq);
    auto topo_end_time = chrono::steady_clock::now();
    auto topo_cost_time = time_cost(topo_start_time, topo_end_time);
    int total_dsu_num = 0, total_topo_unreachable_num = 0;

    while(true)
    {
        int cnt_use_purning = 0;
        ++iterations;
        int dsu_num = 0, topo_unreachable_num = 0;
        vector<vector<int>> vertex_pairs(data_nodes_size);  // sampled vertex pairs
        for (int i=0; i < toSample; ++i)    // counter for reachable pairs
        {
            int seqS = randomNode_temporal(data_nodes_size);
            int seqT = randomNode_temporal(data_nodes_size);
            while(seqS == seqT) {
                seqT = randomNode_temporal(data_nodes_size);
            }
            if(dsu.find(seqS) == dsu.find(seqT)) {
                dsu_num++;
                if(topo_seq[seqS].first <= (topo_seq[seqT].second==INT32_MAX ? topo_seq[seqT].first : topo_seq[seqT].second)) {
                    vertex_pairs[seqS].push_back(seqT);
                } else {
                    topo_unreachable_num++;
                }
            }
        }
        total_dsu_num += dsu_num;
        total_topo_unreachable_num += topo_unreachable_num;
        int sample_seqS_cnt = 0;
        for(auto& vp : vertex_pairs)
        {
            sample_seqS_cnt += vp.size()>0;
        }
        auto calculation_start_time = chrono::steady_clock::now();
        #pragma omp parallel for
        for(int seqS=0; seqS<data_nodes_size; ++seqS)
        {
            if(vertex_pairs[seqS].size() == 0) continue;
            int cur_thread_num = omp_get_thread_num();
            int max_dis = 0, sampled_nodes_num = vertex_pairs[seqS].size(), visited_sample_nodes = 0;
            unordered_map<int, bool> visited_sampled;
            for(auto tmp : vertex_pairs[seqS])
            {
                visited_sampled.emplace(tmp, false);
            }
            Cmp_Node& cmp_root_node = cmp_nodes[seqS];
            vector<int>& distance = t_distance[cur_thread_num];
			vector<int>& cmp_distance = t_cmp_distance[cur_thread_num];
			vector<float>& tmp_sigma = t_tmp_sigma[cur_thread_num];
			vector<float>& total_delta = t_total_delta[cur_thread_num];
			vector<vector<int>>& pre = t_pre[cur_thread_num];
            
            fill(distance.begin(), distance.end(), -1);
            fill(cmp_distance.begin(), cmp_distance.end(), -1);
            fill(tmp_sigma.begin(), tmp_sigma.end(), 0);
            fill(total_delta.begin(), total_delta.end(), 0);
            fill(pre.begin(), pre.end(), vector<int>());
            
            unordered_map<int, int> &total_sigma = tt_total_sigma[cur_thread_num];
            unordered_set<int> &visited_cmp_node = t_visited_cmp_node[cur_thread_num];
            unordered_map<int, int> &cmp_sigma = t_cmp_sigma[cur_thread_num];
            
            unordered_map<int, unordered_map<int, float>> &node_contribs = t_node_contribs[cur_thread_num];
			deque<int> q;
			stack<int> stk;
            total_sigma.clear();
            visited_cmp_node.clear();
            cmp_sigma.clear();
            node_contribs.clear();
            
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
                if( visited_sample_nodes==sampled_nodes_num && cmp_distance[u]>max_dis ) {
                    break;
                }
				for (int j = 0; j < cmp_nodes[u].adj_edge.size(); j++)
				{
					int u_ = cmp_nodes[u].adj_edge[j].destination_seq;
                    int u__seq = -1;
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
                            u__seq = cmp_nodes[u_].identity;
							if (distance[u__seq] == -1)
							{
								distance[u__seq] = cmp_distance[u_];
							}
							if (!visited_cmp_node.count(u_))
							{
								q.push_back(u_);
								visited_cmp_node.insert(u_);
							}
						}
						pre[u_].push_back(u);
					}
                    if(u__seq!=-1 && visited_sampled.count(u__seq)) {
                        if(!visited_sampled[u__seq]) {
                            cnt_use_purning++;
                            visited_sampled[u__seq] = true;
                            visited_sample_nodes++;
                            if(cmp_distance[u_] > max_dis) {
                                max_dis = cmp_distance[u_];
                            }
                        }
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
            float dis;
            while (!stk.empty())
			{
				int u_ = stk.top();
				stk.pop();
				for (int j = 0; j < pre[u_].size(); j++)
				{
					if (pre[u_][j] == cmp_root_node.seq) continue;
					pair<int, int> ut = {pre[u_][j] , cmp_nodes[pre[u_][j]].timestamps.front()};
					if (cmp_nodes[ut.first].identity == cmp_nodes[u_].identity)
					{
                        dis = cmp_nodes[ut.first].timestamps.size() * 1.0 / cmp_nodes[u_].timestamps.size();
						for(auto &id_tbc : node_contribs[u_]) {
                            node_contribs[ut.first][id_tbc.first] += id_tbc.second * dis;
                        }
					}
					else
					{
                        int seq_tmp = cmp_nodes[u_].identity;
						if (visited_sampled.count(seq_tmp) && cmp_distance[u_] == distance[seq_tmp])
						{
                            node_contribs[ut.first][u_] += ((float)(cmp_sigma[ut.first] * cmp_nodes[u_].timestamps.size() * cmp_nodes[ut.first].timestamps.size()) / total_sigma[cmp_nodes[u_].identity]);
                        }
						dis = cmp_sigma[ut.first] * 1.0 * cmp_nodes[ut.first].timestamps.size() / cmp_sigma[u_];
                        for(auto &id_tbc : node_contribs[u_]) {
                            node_contribs[ut.first][id_tbc.first] += id_tbc.second * dis;

                        }
					}
				}
			}
            
            unordered_map<int, unordered_map<int, double>> importance2;
            for(auto &item1 : node_contribs) {
                for(auto &item2 : item1.second) {
                    if(item1.first == item2.first) {
                        continue;
                    }
                    importance2[cmp_nodes[item1.first].identity][cmp_nodes[item2.first].identity] += item2.second;
                }
            }

            #pragma omp critical
            {
                for(auto &nt_cb : importance2)
                {
                    for(auto &cb : nt_cb.second) {
                        sumContribs[nt_cb.first] += cb.second;
                        sumContribsSquared[nt_cb.first] += cb.second * cb.second;
                    }
                }
            }
        }

        auto calculation_end_time = chrono::steady_clock::now();
        auto calculation_cost_time = time_cost(calculation_start_time, calculation_end_time);
        std::cout << "calculation time = " << calculation_cost_time << endl;
        total_calculation_cost_time += calculation_cost_time;

        auto rademacher_start_time = chrono::steady_clock::now();

        const double exponentDenominator = 2 * (sampleSize * sampleSize);
        // cout << "exponentDenominator = " << exponentDenominator << endl;
        RunningObjFuncSumExponents runningObjFuncSumExponents;
        // cout << "contribs:\n";
        for(int cur=0; cur<data_nodes_size; ++cur)
        {
            scoreData[cur] = sumContribs[cur];
            int hash = 0;
            double objFuncSumExponent = sumContribsSquared[cur];
            hash ^= cur + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            runningObjFuncSumExponents += make_pair(hash, objFuncSumExponent / exponentDenominator);
        }
        vector<double> objFuncSumExponents = runningObjFuncSumExponents.getVector();
        supDeviationBound = NetworKit::getSupDeviationBound2(setting.delta, sampleSize, setting.penaltyFactor, objFuncSumExponents, minimizationResults);
        double epsilon = setting.epsilon;
        if(supDeviationBound > epsilon)
        {
            // compute new sample size and iterate
            toSample = sampleSize;
            assert(minimizationResults.first < epsilon);
            const double mylog = log(2 / setting.delta);
            const double epsilon_square = epsilon * epsilon;
            const double radebound_square = minimizationResults.first * minimizationResults.first;
            const double twentyseven_radebound_square = 27 * radebound_square;
            const double next_multiplier = mylog / (6 * (epsilon_square - 2 * epsilon * minimizationResults.first + radebound_square));
            const double next_1st_term = 2 + 8 * epsilon;
            const double next_2nd_term_multiplier = sqrt(48 * minimizationResults.first + 1 + 8 * epsilon + 16 * epsilon_square);
            const double atan2_1st_arg = 12 * sqrt(3) * fabs(-1 + 2 * (minimizationResults.first + epsilon)) * sqrt(- (twentyseven_radebound_square - epsilon_square * (1 + 16 * epsilon) - minimizationResults.first * (1 + 18 * epsilon)));
            const double atan2_2nd_arg = - ( -1 - 12 * epsilon + 8 * ( twentyseven_radebound_square + (21 - 8 * epsilon) * epsilon_square + 18 * minimizationResults.first * (1 + epsilon)));
            const double atan2_res = atan2(atan2_1st_arg, atan2_2nd_arg);
            const double theta = atan2_res / 3;
            const double cos_theta = cos(theta);
            const double sin_theta = sin(theta);
            const array<uint64_t, 3> roots = {
                (uint64_t) ceil(next_multiplier * (next_1st_term - 2 * next_2nd_term_multiplier * cos_theta)),
                (uint64_t) ceil(next_multiplier * (next_1st_term + next_2nd_term_multiplier * (cos_theta + sqrt(3) * sin_theta))),
                (uint64_t) ceil(next_multiplier * (next_1st_term + next_2nd_term_multiplier * (cos_theta - sqrt(3) * sin_theta)))};
            sampleSize = *(max_element(roots.begin(), roots.end()));
            for (auto root : roots) {
                if (root <= toSample || root == sampleSize) {
                    continue;
                }
                const double alpha = mylog / (mylog + sqrt(mylog * (2 * root * minimizationResults.first + mylog)));
                const double first_term = sqrt(mylog / (2 * root));
                const double second_term = minimizationResults.first / (1 - alpha);
                const double third_term = mylog / (2 * root * alpha * (1 - alpha));
                const double next_deviation = first_term + second_term + third_term;
                if (next_deviation <= epsilon && root < sampleSize) {  // if not only one root satisified the requirement, choose the smallest one
                    sampleSize = root;
                }
            }
            assert(sampleSize > toSample);
            toSample = sampleSize - toSample;
            auto rademacher_end_time = chrono::steady_clock::now();
            auto rademacher_cost_time = time_cost(rademacher_start_time, rademacher_end_time);
            total_rademacher_cost_time += rademacher_cost_time;
        } else {
            auto rademacher_end_time = chrono::steady_clock::now();
            auto rademacher_cost_time = time_cost(rademacher_start_time, rademacher_end_time);
            total_rademacher_cost_time += rademacher_cost_time;
            break;
        }
    }
    std::cout << "total time = " << total_calculation_cost_time + total_rademacher_cost_time
              << "\ntotal_calculation_cost_time = " << total_calculation_cost_time
              << "\ntotal_rademacher_cost_time = " << total_rademacher_cost_time
              << "\nsampleSize = "<< sampleSize << endl;
            //   << "\ntotal_unreachable_pairs = " << total_unreachable_pairs << endl;
    
    ofstream out("./bc-results/"+dataName+"-"+to_string(setting.epsilon).substr(0, 5)+".csv");
    if(out.is_open()) {
        for(int i=0; i<data_nodes_size; ++i)
        {
            bcResult[i] = scoreData[i] / sampleSize;
            out << i << "," << seq_to_id[i] << "," << bcResult[i] << "," << scoreData[i] << endl;
        }
        out.close();
    } else {
        clog << "write bc results failed!\n";
    }
    return bcResult;
}
