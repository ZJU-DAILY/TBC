#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <malloc.h>
#include <sys/types.h>
#include <dirent.h>
#include <chrono>
#include <iomanip>
#include "Transformer.h"
#include "Algorithm.h"

using namespace std;

int main()  
{
    Settings setting;
    vector<double> epsilons({0.03, 0.025, 0.02, 0.015, 0.01, 0.005});  //0.03, 0.025, 0.02, 0.015, 0.01, 0.005
    std::vector<std::string> files;
    string relativePath = "../../data-set/";
    // DIR* dir = opendir(relativePath.c_str());
    // if (dir == NULL)
    // {
    //     printf("dir == NULL");
    // }
    // struct dirent* entry;
    // while ( (entry=readdir(dir)) != NULL)
    // {
    //     std::string str(entry->d_name);
    //     if(str.size() < 4)
    //         continue;
    //     files.push_back(str);
    // }
    // closedir(dir);

    // files.push_back("edit-sewiki.txt");
    // files.push_back("edit-itwikivoyage.txt");
    files.push_back("infectious.txt");
    // files.push_back("sx-mathoverflow.txt");
    // files.push_back("superuser.txt");
    // files.push_back("wiki-talk-temporal.txt");
    for(auto file : files)
    {
        string filePath = relativePath + file;
        string dataName = file.substr(0, file.size()-4);
        cout << "START " << dataName << "*****************************************" << endl;
        Transformer transformer = Transformer(filePath);

        auto &nodes = transformer.data_graph.graph_data;
        auto &edge_set = transformer.data_graph.edge_set;
        Dsu dsu(nodes.size());
        for(auto &e : edge_set) {
            dsu.unite(e.first.first, e.first.second);
        }
        transformer.doTransform(setting.isStrict);
        if(setting.needCompress)
        {
            transformer.doCompress();
            for(auto ep : epsilons) {
                time_t now = std::chrono::system_clock::to_time_t(chrono::system_clock::now());
                clog << put_time(localtime(&now), "%F %T") <<  endl;
                cout << "\nsample algorithm start with epsilon = " << ep << endl;
                clog << dataName << " + " << ep << "start\n";
                setting.epsilon = ep;
                sample_compress_v_v_parallel_purning_experiment(transformer, setting, dataName, dsu);
            }
            cout << "*****************************************END " << dataName << "\n\n";
            malloc_trim(0);
            clog << dataName << "  finished\n";
        }
    }
}