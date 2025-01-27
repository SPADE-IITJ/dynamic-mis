// #include <omp.h>
// #include <iostream>
// #include <vector>
// #include <fstream>
// #include <sstream>
// #include <unordered_map>
// #include <unordered_set>
// #include <chrono>
// #include <algorithm>
// #include <filesystem>
// #include <iomanip>
// #include <math.h>
// #include <set>
// #include <stdbool.h>
// #include <stdio.h>
// // #define num_threads 64
// int n_threads = 75;
// using namespace std;
// vector<double> time_per_thread;
// double maxi_time = 0;
// // using namespace fs = std::filesystem;
// // vector<int> mis;
// class Graph
// {
// public:
//     std::vector<int> edges;        // stores destination nodes of edges
//     std::vector<long int> indices; // stores index in 'edges' for the start of each node's edge list
//     std::vector<int> degree;       // stores the degree of each node
//     std::vector<bool> nodeStatus;
//     std::vector<vector<int>> tvbl;
//     std::vector<int> csredges;
//     std::vector<int> csrindices;
    
//     int maxNode = -1;

//     vector<bool> isolated;
//     bool isIndependentSet(vector<int> mis)
//     {
//         for (int v : mis)
//         {
//             for (int u : mis)
//             {
//                 if (v != u && this->isAdjacent(v, u))
//                 {
//                     return false; // Not an independent set
//                 }
//             }
//         }
//         return true;
//     }
//     // This method is used to build the CSR from the edge list
//     void buildCSR(const std::vector<std::pair<int, int>> &edgeList)
//     {
//         // Find the maximum node value
//         for (const auto &edge : edgeList)
//         {
//             maxNode = std::max({maxNode, edge.first, edge.second});
//         }

//         // Resize the vectors to accommodate the maximum node value
//         indices.assign(maxNode + 1, 0); // We need maxNode + 2 to include 0 and maxNode indices
//         edges.resize(edgeList.size() * 2, 0);
//         degree.assign(maxNode + 1, 0); // Since it's an undirected graph, each edge will be added twice
//         nodeStatus.assign(maxNode + 1, false);
//         isolated.assign(maxNode + 1, false);

//         // Calculate the degree of each node
//         for (const auto &edge : edgeList)
//         {
//             degree[edge.first]++;
//             degree[edge.second]++;
//         }

//         // Populate the indices vector
//         for (int i = 1; i <= maxNode; ++i)
//         {
//             indices[i] = indices[i - 1] + degree[i - 1];
//         }
//         std::vector<long int> currInd = indices;

//         // Populate the edges vector
//         std::vector<int> indexCount = degree; // Temporary array to keep track of the current index position for each node
//         for (const auto &edge : edgeList)
//         {
//             edges[currInd[edge.first]++] = edge.second;
//             edges[currInd[edge.second]++] = edge.first;
//         }
//         for (int i = 0; i < degree.size(); i++)
//         {
//             if (degree[i] == 0)
//             {
//                 nodeStatus[i] = true;
//                 isolated[i] = true;
//             }
//         }

//         cout << "CSR made" << endl;
//         cout << "Degree->" << degree.size() << endl;
//         cout << "Indices->" << indices.size() << endl;
//         cout << "Edges->" << edges.size() / 2 << endl;
//     }

//     // Checks if a node 'u' is adjacent to node 'v'
//     bool isAdjacent(int u, int v) const
//     {
//         if (u >= csrindices.size() - 1 || v >= csrindices.size() - 1)
//             return false;

//         for (int i = csrindices[u]; i < csrindices[u + 1]; ++i)
//         {
//             if (csredges[i] == v)
//                 return true;
//         }
//         return false;
//     }

//     // Get the neighbors of a node
//     std::vector<int> getNeighbors(int u) const
//     {
//         std::vector<int> neighbors;
//         if (u < indices.size() - 1)
//         {
//             for (int i = csrindices[u]; i < csrindices[u + 1]; ++i)
//             {
//                 neighbors.push_back(csredges[i]);
//             }
//         }
//         if (u == indices.size() - 1)
//         {
//             for (int i = csrindices[u]; i < csredges.size(); i++)
//             {
//                 neighbors.push_back(csredges[i]);
//             }
//         }
//         return neighbors;
//     }

//     std::vector<int> getMNeighbors(int u) const
//     {
//         std::vector<int> neighbors;
//         u--;
//         if (u < indices.size() - 1)
//         {
//             for (int i = indices[u]; i < indices[u + 1]; ++i)
//             {
//                 if (edges[i] != -1)
//                 {

//                     neighbors.push_back(edges[i]);
//                 }
                
//             }
//         }
//         if (u == indices.size() - 1)
//         {
//             for (int i = indices[u]; i < edges.size(); i++)
//             {
//                 // neighbors.push_back(edges[i]);
//                 if (edges[i] != -1)
//                 {

//                     neighbors.push_back(edges[i]);
//                 }
//             }
//         }
//         return neighbors;
//     }
// };

// void PDMI_INC(Graph &graph, vector<int> &mis, pair<int, int> edge)
// {
//     bool firstIn = false, secondIn = false;
//     int it = -1;
//     bool firstInLocal, secondInLocal;
//     double beg, en;
//     time_per_thread.clear();
//     time_per_thread.resize(n_threads);
// #pragma omp parallel num_threads(n_threads) shared(firstIn, secondIn, it) private(firstInLocal, secondInLocal)
//     {
//         beg = omp_get_wtime();
//         int tid = omp_get_thread_num();
//         int start = tid * (mis.size() / n_threads);
//         int end = (tid == n_threads - 1) ? mis.size() : (tid + 1) * (mis.size() / n_threads);

//         firstInLocal = false, secondInLocal = false;

//         for (int i = start; i < end; i++)
//         {
//             if (edge.first == mis[i])
//             {
//                 firstInLocal = true;
//                 it = i;
//             }
//             if (edge.second == mis[i])
//             {
//                 secondInLocal = true;
//             }
//         }

// #pragma omp atomic update
//         firstIn |= firstInLocal;
// #pragma omp atomic update
//         secondIn |= secondInLocal;

//         en = omp_get_wtime();
//         time_per_thread[tid] = (en - beg);
//     }
//     if (firstIn && secondIn && it != -1)
//     {

//         mis.erase(mis.begin() + it);
//         vector<int> neighbours = graph.getNeighbors(edge.first);
//         int p = 0;
//         for (int k = 0; k < neighbours.size(); k++)
//         {
//             bool shouldWeAdd = true;
//             // printf("ji");
//             for (int l = 0; l < graph.getNeighbors(neighbours[k]).size(); l++)
//             {
//                     bool shouldweAddLocal = true;
// #pragma omp parallel num_threads(n_threads) shared(graph, neighbours, mis, shouldWeAdd) private(shouldweAddLocal)
//                 {
//                     beg = omp_get_wtime();
//                     int tid = omp_get_thread_num();
//                     int start = tid * (mis.size() / n_threads);
//                     int end = (tid == n_threads - 1) ? mis.size() : (tid + 1) * (mis.size() / n_threads);
//                     for (int i = start; i < end; i++)
//                     {
//                         if (graph.getNeighbors(neighbours[k])[l] == mis[i])
//                         {
//                             shouldweAddLocal = false;
                           
//                         }
//                     }
                    
// #pragma omp atomic update
//         shouldWeAdd |= shouldweAddLocal;
//         end = omp_get_wtime();
//         time_per_thread[tid] += end - beg;
//                 }
//                 if (shouldWeAdd)
//                 {
//                     mis.push_back(neighbours[k]);
//                 }
//             }
//         }
//     }
// }

// int main(int argc, char *argv[])
// {
//     if (argc < 7)
//     {
//         std::cerr << "Usage: " << argv[0] << " <CSR> <MIS> <Deletion_folder> <Insertion_folder> <num_threads> <Num_batch for deletion> <Num_batch for insertion>" << std::endl;
//         return 1;
//     }
//     cout << "----------------------------------------Start-----------------------------------------" << endl;
    
//     cout << "----------------------------------------INC_BAT------------------------------------------" << endl;
//     std::string filename1 = argv[1];
//     std::string filename2 = argv[2];
//     Graph graph;

//     graph = readCSR(filename1, filename2);
//     std::vector<int> maximalIndependentSet;
//     std::string line;
//     // Read MIS------------------------------
//     vector<bool> nodeStats(graph.maxNode + 1, false);

//     std::ifstream inputFile(argv[2]);
//     if (inputFile.is_open())
//     {
//         while (std::getline(inputFile, line))
//         {
//             maximalIndependentSet.push_back(std::stoi(line) + 1);
//             nodeStats[std::stoi(line) + 1] = true;
//         }
//         inputFile.close();
//     }
//     else
//     {
//         std::cerr << "Unable to open the file for reading" << std::endl;
//     }

//     vector<int> misforInsertionINC;
//     std::copy(maximalIndependentSet.begin(), maximalIndependentSet.end(), std::back_inserter(misforInsertionINC));
//     vector<int> misforInsertionBAT;
//     std::copy(maximalIndependentSet.begin(), maximalIndependentSet.end(), std::back_inserter(misforInsertionBAT));


//     vector<bool> nodeStatsforInsertion(nodeStats.size(), false);
//     for (int i = 0; i < nodeStats.size(); i++)
//     {
//         if (nodeStats[i])
//         {
//             nodeStatsforInsertion[i] = true;
//         }
//     }

//     std::string folderPath(argv[3]);
//     std::vector<std::pair<int, int>> edgeList;
//     std::vector<std::vector<std::pair<int, int>>> edgeListsForDeletion;
//     for (const auto &entry : std::filesystem::directory_iterator(folderPath))
//     {
//         if (entry.is_regular_file())
//         {
//             std::ifstream file(entry.path());
//             std::string line;
//             while (std::getline(file, line))
//             {
//                 std::istringstream iss(line);
//                 int src, dest;
//                 if (iss >> src)
//                 {
//                     if (!(iss >> dest))
//                     {
//                         dest = -1;
//                     }
//                     edgeList.push_back(std::make_pair(src, dest));
//                 }
//             }
//             edgeListsForDeletion.push_back(edgeList);
//             edgeList.clear();
//         }
//     }

//     std::string folderPath1(argv[4]);
//     std::vector<std::vector<std::pair<int, int>>> edgeListsForInsertion;
//     for (const auto &entry : std::filesystem::directory_iterator(folderPath1))
//     {
//         if (entry.is_regular_file())
//         {
//             std::ifstream file(entry.path());
//             std::string line;
//             while (std::getline(file, line))
//             {
//                 std::istringstream iss(line);
//                 int src, dest;
//                 if (iss >> src)
//                 {
//                     if (!(iss >> dest))
//                     {
//                         dest = -1;
//                     }
//                     edgeList.push_back(std::make_pair(src, dest));
//                 }
//             }
//             edgeListsForInsertion.push_back(edgeList);
//             edgeList.clear();
//         }
//     }
//     int initialCard = maximalIndependentSet.size();
//     n_threads = stoi(argv[5]);
//     double startTime = 0, endTime = 0, time_cumm = 0;
//     int num_deletion_batches = std::stoi(argv[6]);
//     int num_insertion_batches = std::stoi(argv[7]);
//     cout<< "Cores USed-: "<< n_threads<<endl;
//     cout << "Graph Name->" << filename1 << endl;
//     cout <<"Batch Size->"<<edgeListsForInsertion[0].size() <<" | Number of batches->"<<num_insertion_batches<< " | Incremental Edges->" << edgeListsForInsertion[0].size() * num_insertion_batches << endl;
//     cout <<"Batch Size->"<<edgeListsForDeletion[0].size() <<" | Number of batches->"<<num_deletion_batches<< " | Deletion Edges->" << edgeListsForInsertion[0].size() * num_insertion_batches << endl;
//     cout << "MIS set size->" << maximalIndependentSet.size() << endl;


// //----------------------------------------------INCREMENTAL---------------------------------------------------------

//     auto start = omp_get_wtime();
//     double time_cum = 0;
//     time_per_thread.clear();
//     time_per_thread.resize(n_threads, 0.0);

//     time_cum = 0;
//     start = omp_get_wtime();
//     for (int i = 0; i < num_insertion_batches; i++)
//     {
//         for (int j = 0; j < edgeListsForInsertion[i].size(); j++)
//         {
//             cout<<edgeListsForInsertion[i][j].first<<endl;
//             PDMI_INC(graph, misforInsertionINC, edgeListsForInsertion[i][j]);
//             double maxi = 0;
//             for (int k = 0; k < time_per_thread.size(); k++)
//             {
//                 // calculate the max element
//                 if (maxi < time_per_thread[k])
//                 {
//                     maxi = time_per_thread[k];
//                 }
//             }
//             time_cumm += maxi;
//         }
//     }
//     auto end = omp_get_wtime();
//     cout << endl;
//     cout << "-------------------------------INCREMENTAL------------------------------" << endl;
//     cout << "Cardinality changed by->(Insertion)" << abs(initialCard - static_cast<int>(misforInsertionINC.size())) << endl;
//     cout << "Total Incremental Edges Insertion time: " << fixed << setprecision(5) << time_cumm * 1000 << " milliseconds" << endl;
//     cout << "Average Incremental Edges Insertion time: " << fixed << setprecision(5) << time_cumm / num_deletion_batches * 1000 << " milliseconds" << endl;



// //----------------------------------------------BATCH_DYNAMIC---------------------------------------------------------



//     time_per_thread.clear();
//     time_per_thread.resize(n_threads, 0.0);

//         time_cumm = 0;
//         for (int i = 0; i < num_insertion_batches; i++)
//         {
//             int block, tid, start, end, it;
//             omp_set_num_threads(n_threads);
//             block = ceil(static_cast<float>(edgeListsForInsertion[i].size()) / n_threads);
//     #pragma omp parallel private(tid, start, end, it) shared(misforInsertionBAT)
//             {
//                 startTime = omp_get_wtime();
//                 tid = omp_get_thread_num();
//                 start = tid * block;
//                 end = std::min((tid + 1) * block, static_cast<int>(edgeListsForInsertion[i].size()));

//                 for (int j = start; j < end; j++)
//                 {
//                     bool firstIn = false, secondIn = false;
//                     it = 0;

//                     for (int k = 0; k < misforInsertionBAT.size(); k++)
//                     {
//                         if (edgeListsForInsertion[i][j].first == misforInsertionBAT[k])
//                         {
//                             {
//                                 firstIn = true;
//                                 it = k;
//                             }
//                         }
//                         if (edgeListsForInsertion[i][j].second == misforInsertionBAT[k])
//                         {
//                             {
//                                 secondIn = true;
//                             }
//                         }
//                     }

//                     if (firstIn && secondIn)
//                     {
//                         misforInsertionBAT[it] = -1;
//                                vector<int> neighbours=graph.getNeighbors(edgeListsForInsertion[i][j].first);
//             int p=0;
//                     for(int k=0;k<neighbours.size();k++)
//                     {
//                     bool shouldWeAdd=true;
//                         for(int l=0;l<graph.getNeighbors(neighbours[k]).size();l++){
//                             for(int m=0;m<misforInsertionBAT.size();m++){
//                                 if(misforInsertionBAT[m]==graph.getNeighbors(neighbours[k])[l]){
//                                     shouldWeAdd=false;
//                                     break;
//                                 }
//                             }
//                         }
//                     if(shouldWeAdd){
//                         misforInsertionBAT.push_back(neighbours[k]);
//                     }
//                     }
//                     }
//                 }
//                 endTime = omp_get_wtime();
//                 time_per_thread[tid] = endTime - startTime;
//             }
//             double maxi = 0;
//             for (int k = 0; k < time_per_thread.size(); k++)
//             {
//                 // calculate the max element
//                 if (maxi < time_per_thread[k])
//                 {
//                     maxi = time_per_thread[k];
//                 }
//             }
//             time_cumm += maxi;
//         }
//         int ca = 0;
//         for (int i = 0; i < misforInsertionBAT.size(); i++)
//         {
//             if (misforInsertionBAT[i] == -1)
//             {
//                 ca++;
//             }
//         }
//         cout << endl;
//         cout << "-------------------------------BATCHDYNAMIC------------------------------" << endl;
//         cout << "Batch size of one small batch is-> " << edgeListsForInsertion[0].size() << endl;
//         cout << "Change in Cardinality " << abs(ca) << endl;
//         cout << "Parallel Edges Insertion time: " << fixed << setprecision(2) << time_cumm * 1000 << " milliseconds" << endl;
//         cout << "Average Parallel Edges Insertion time: " << fixed << setprecision(2) << time_cumm / num_insertion_batches * 1000 << " milliseconds" << endl;

//     cout << "---------------------------------------END-----------------------------------------" << endl<<endl<<endl;
// }



/**
 * @file inc_Bat.cpp
 * @brief Implementation PDMI_INC and PDMI_BAT algorithms
 */

#include <omp.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <math.h>
#include <set>
#include <stdbool.h>
#include <stdio.h>

// Global variables
int n_threads = 75;
using namespace std;
vector<double> time_per_thread;

/**
 * @class Graph
 * @brief Represents a graph using Compressed Sparse Row (CSR) format
 */
class Graph
{
public:
    std::vector<int> csredges;
    std::vector<int> csrindices;
    
    int maxNode = -1;


    /**
     * @brief Checks if the given set of nodes forms an independent set
     * @param mis Vector of nodes to check
     * @return true if the set is independent, false otherwise
     */
    bool isIndependentSet(vector<int> mis)
    {
        for (int v : mis)
        {
            for (int u : mis)
            {
                if (v != u && this->isAdjacent(v, u))
                {
                    return false; // Not an independent set
                }
            }
        }
        return true;
    }

    /**
     * @brief Checks if two nodes are adjacent in the graph
     * @param u First node
     * @param v Second node
     * @return true if the nodes are adjacent, false otherwise
     */
    bool isAdjacent(int u, int v) const
    {
        if (u >= csrindices.size() - 1 || v >= csrindices.size() - 1)
            return false;

        for (int i = csrindices[u]; i < csrindices[u + 1]; ++i)
        {
            if (csredges[i] == v)
                return true;
        }
        return false;
    }

    /**
     * @brief Gets the neighbors of a given node
     * @param u Node to get neighbors for
     * @return Vector of neighboring nodes
     */
    std::vector<int> getNeighbors(int u) const
    {
        std::vector<int> neighbors;
        if (u < csrindices.size() - 1)
        {
            for (int i = csrindices[u]; i < csrindices[u + 1]; ++i)
            {
                neighbors.push_back(csredges[i]);
            }
        }
        if (u == csrindices.size() - 1)
        {
            for (int i = csrindices[u]; i < csredges.size(); i++)
            {
                neighbors.push_back(csredges[i]);
            }
        }
        return neighbors;
    }
};


Graph readCSR(string filename)

{

    Graph graph;
    std::vector<int> indices, edges;
    std::ifstream file(filename);
    if (file.is_open())
    {
        int nodes, edgesCount;
        file >> nodes >> edgesCount;
        std::string line;
        bool x = 1;
        // Read the indices array
        while (std::getline(file, line))
        {
            if (line == "indices->")
            {
                while (std::getline(file, line))
                {
                    if (line == "list->")
                    {
                        x = 0;
                        break;
                    }
                    indices.push_back(std::stoi(line));
                }
                while (std::getline(file, line))
                {
                    edges.push_back(std::stoi(line) + 1);
                }
                break;
            }
        }
    }
    graph.maxNode = indices.size();
    graph.csrindices.resize(indices.size() + 1, 0);
    graph.csredges.resize(edges.size(),0);
    graph.csrindices[0] = 0;
    for (int i = 1; i <= indices.size(); i++)
    {
        graph.csrindices[i] = indices[i - 1];
    }

    for (int i = 0; i < edges.size(); i++)
    {
        graph.csredges[i] = edges[i] ;
    }

    cout << "Node->" << graph.maxNode + 1 << endl;
    cout << "Edges->" << graph.csredges.size() << endl;

    return graph;
}
/**
 * @brief Performs incremental update of Maximal Independent Set
 * @param graph The graph structure
 * @param mis Vector representing the Maximal Independent Set
 * @param edge Pair of nodes representing the new edge
 */
void PDMI_INC(Graph &graph, vector<int> &mis, pair<int, int> edge)
{
    bool firstIn = false, secondIn = false;
    int it = -1;
    bool firstInLocal, secondInLocal;
    double beg, en;
    time_per_thread.clear();
    time_per_thread.resize(n_threads);

    // Parallel section to check if nodes are in MIS
    #pragma omp parallel num_threads(n_threads) shared(firstIn, secondIn, it) private(firstInLocal, secondInLocal)
    {
        beg = omp_get_wtime();
        int tid = omp_get_thread_num();
        int start = tid * (mis.size() / n_threads);
        int end = (tid == n_threads - 1) ? mis.size() : (tid + 1) * (mis.size() / n_threads);

        firstInLocal = false, secondInLocal = false;

        for (int i = start; i < end; i++)
        {
            if (edge.first == mis[i])
            {
                firstInLocal = true;
                it = i;
            }
            if (edge.second == mis[i])
            {
                secondInLocal = true;
            }
        }

        #pragma omp atomic update
        firstIn |= firstInLocal;
        #pragma omp atomic update
        secondIn |= secondInLocal;

        en = omp_get_wtime();
        time_per_thread[tid] = (en - beg);
    }

    // Update MIS if both nodes of the new edge are in MIS
    if (firstIn && secondIn && it != -1)
    {
        mis.erase(mis.begin() + it);
        vector<int> neighbours = graph.getNeighbors(edge.first);
        int p = 0;
        for (int k = 0; k < neighbours.size(); k++)
        {
            bool shouldWeAdd = true;
            for (int l = 0; l < graph.getNeighbors(neighbours[k]).size(); l++)
            {
                bool shouldweAddLocal = true;
                #pragma omp parallel num_threads(n_threads) shared(graph, neighbours, mis, shouldWeAdd) private(shouldweAddLocal)
                {
                    beg = omp_get_wtime();
                    int tid = omp_get_thread_num();
                    int start = tid * (mis.size() / n_threads);
                    int end = (tid == n_threads - 1) ? mis.size() : (tid + 1) * (mis.size() / n_threads);
                    for (int i = start; i < end; i++)
                    {
                        if (graph.getNeighbors(neighbours[k])[l] == mis[i])
                        {
                            shouldweAddLocal = false;
                        }
                    }
                    
                    #pragma omp atomic update
                    shouldWeAdd |= shouldweAddLocal;
                    end = omp_get_wtime();
                    time_per_thread[tid] += end - beg;
                }
                if (shouldWeAdd)
                {
                    mis.push_back(neighbours[k]);
                }
            }
        }
    }
}

/**
 * @brief Main function to process graph and perform MIS operations
 * @param argc Number of command-line arguments
 * @param argv Array of command-line arguments
 * @return 0 on successful execution
 */
int main(int argc, char *argv[])
{
    // Check for correct number of command-line arguments
    if (argc < 7)
    {
        std::cerr << "Usage: " << argv[0] << " <CSR> <MIS> <Deletion_folder> <Insertion_folder> <num_threads> <Num_batch for deletion> <Num_batch for insertion>" << std::endl;
        return 1;
    }

    cout << "----------------------------------------Start-----------------------------------------" << endl;
    cout << "----------------------------------------INC_BAT------------------------------------------" << endl;

    // Read input files and initialize graph
    std::string filename1 = argv[1];
    std::string filename2 = argv[2];
    Graph graph;

    graph = readCSR(filename1);
    std::vector<int> maximalIndependentSet;
    std::string line;

    // Read Maximal Independent Set

    std::ifstream inputFile(argv[2]);
    if (inputFile.is_open())
    {
        while (std::getline(inputFile, line))
        {
            maximalIndependentSet.push_back(std::stoi(line) + 1);
        }
        inputFile.close();
    }
    else
    {
        std::cerr << "Unable to open the file for reading" << std::endl;
    }

    // Initialize MIS for incremental and batch processing
    vector<int> misforInsertionINC(maximalIndependentSet);
    vector<int> misforInsertionBAT(maximalIndependentSet);

    // Read edge lists for deletion and insertion
    std::string folderPath(argv[3]);
    std::vector<std::pair<int, int>> edgeList;
    std::vector<std::vector<std::pair<int, int>>> edgeListsForDeletion;
    for (const auto &entry : std::filesystem::directory_iterator(folderPath))
    {
        if (entry.is_regular_file())
        {
            std::ifstream file(entry.path());
            std::string line;
            while (std::getline(file, line))
            {
                std::istringstream iss(line);
                int src, dest;
                if (iss >> src)
                {
                    if (!(iss >> dest))
                    {
                        dest = -1;
                    }
                    edgeList.push_back(std::make_pair(src, dest));
                }
            }
            edgeListsForDeletion.push_back(edgeList);
            edgeList.clear();
        }
    }

    std::string folderPath1(argv[4]);
    std::vector<std::vector<std::pair<int, int>>> edgeListsForInsertion;
    for (const auto &entry : std::filesystem::directory_iterator(folderPath1))
    {
        if (entry.is_regular_file())
        {
            std::ifstream file(entry.path());
            std::string line;
            while (std::getline(file, line))
            {
                std::istringstream iss(line);
                int src, dest;
                if (iss >> src)
                {
                    if (!(iss >> dest))
                    {
                        dest = -1;
                    }
                    edgeList.push_back(std::make_pair(src, dest));
                }
            }
            edgeListsForInsertion.push_back(edgeList);
            edgeList.clear();
        }
    }

    // Initialize parameters
    int initialCard = maximalIndependentSet.size();
    n_threads = stoi(argv[5]);
    double startTime = 0, endTime = 0, time_cumm = 0;
    int num_deletion_batches = std::stoi(argv[6]);
    int num_insertion_batches = std::stoi(argv[7]);

    // Print initial information
    cout << "Cores Used-: " << n_threads << endl;
    cout << "Graph Name->" << filename1 << endl;
    cout << "Batch Size->" << edgeListsForInsertion[0].size() << " | Number of batches->" << num_insertion_batches << " | Incremental Edges->" << edgeListsForInsertion[0].size() * num_insertion_batches << endl;
    cout << "Batch Size->" << edgeListsForDeletion[0].size() << " | Number of batches->" << num_deletion_batches << " | Deletion Edges->" << edgeListsForInsertion[0].size() * num_insertion_batches << endl;
    cout << "MIS set size->" << maximalIndependentSet.size() << endl;

    //----------------------------------------------INCREMENTAL---------------------------------------------------------

    auto start = omp_get_wtime();
    double time_cum = 0;
    time_per_thread.clear();
    time_per_thread.resize(n_threads, 0.0);

    time_cum = 0;
    start = omp_get_wtime();
    for (int i = 0; i < num_insertion_batches; i++)
    {
        for (int j = 0; j < edgeListsForInsertion[i].size(); j++)
        {
            // cout << edgeListsForInsertion[i][j].first << endl;
            PDMI_INC(graph, misforInsertionINC, edgeListsForInsertion[i][j]);
            double maxi = 0;
            for (int k = 0; k < time_per_thread.size(); k++)
            {
                // calculate the max element
                if (maxi < time_per_thread[k])
                {
                    maxi = time_per_thread[k];
                }
            }
            time_cumm += maxi;
        }
    }
    auto end = omp_get_wtime();
    cout << endl;
    cout << "-------------------------------INCREMENTAL------------------------------" << endl;
    cout << "Cardinality changed by->(Insertion)" << abs(initialCard - static_cast<int>(misforInsertionINC.size())) << endl;
    cout << "Total Incremental Edges Insertion time: " << fixed << setprecision(5) << time_cumm * 1000 << " milliseconds" << endl;
    cout << "Average Incremental Edges Insertion time: " << fixed << setprecision(5) << time_cumm / num_deletion_batches * 1000 << " milliseconds" << endl;

    //----------------------------------------------BATCH_DYNAMIC---------------------------------------------------------

    time_per_thread.clear();
    time_per_thread.resize(n_threads, 0.0);

    time_cumm = 0;
    for (int i = 0; i < num_insertion_batches; i++)
    {
        int block, tid, start, end, it;
        omp_set_num_threads(n_threads);
        block = ceil(static_cast<float>(edgeListsForInsertion[i].size()) / n_threads);
#pragma omp parallel private(tid, start, end, it) shared(misforInsertionBAT)
        {
            startTime = omp_get_wtime();
            tid = omp_get_thread_num();
            start = tid * block;
            end = std::min((tid + 1) * block, static_cast<int>(edgeListsForInsertion[i].size()));

            for (int j = start; j < end; j++)
            {
                bool firstIn = false, secondIn = false;
                it = 0;

                for (int k = 0; k < misforInsertionBAT.size(); k++)
                {
                    if (edgeListsForInsertion[i][j].first == misforInsertionBAT[k])
                    {
                        {
                            firstIn = true;
                            it = k;
                        }
                    }
                    if (edgeListsForInsertion[i][j].second == misforInsertionBAT[k])
                    {
                        {
                            secondIn = true;
                        }
                    }
                }

                if (firstIn && secondIn)
                {
                    misforInsertionBAT[it] = -1;
                    vector<int> neighbours = graph.getNeighbors(edgeListsForInsertion[i][j].first);
                    int p = 0;
                    for (int k = 0; k < neighbours.size(); k++)
                    {
                        bool shouldWeAdd = true;
                        for (int l = 0; l < graph.getNeighbors(neighbours[k]).size(); l++) {
                            for (int m = 0; m < misforInsertionBAT.size(); m++) {
                                if (misforInsertionBAT[m] == graph.getNeighbors(neighbours[k])[l]) {
                                    shouldWeAdd = false;
                                    break;
                                }
                            }
                        }
                        if (shouldWeAdd) {
                            misforInsertionBAT.push_back(neighbours[k]);
                        }
                    }
                }
            }
            endTime = omp_get_wtime();
            time_per_thread[tid] = endTime - startTime;
        }
        double maxi = 0;
        for (int k = 0; k < time_per_thread.size(); k++)
        {
            // calculate the max element
            if (maxi < time_per_thread[k])
            {
                maxi = time_per_thread[k];
            }
        }
        time_cumm += maxi;
    }
    int ca = 0;
    for (int i = 0; i < misforInsertionBAT.size(); i++)
    {
        if (misforInsertionBAT[i] == -1)
        {
            ca++;
        }
    }
    cout << endl;
    cout << "-------------------------------BATCHDYNAMIC------------------------------" << endl;
    cout << "Batch size of one small batch is-> " << edgeListsForInsertion[0].size() << endl;
    cout << "Change in Cardinality " << abs(ca) << endl;
    cout << "Parallel Edges Insertion time: " << fixed << setprecision(2) << time_cumm * 1000 << " milliseconds" << endl;
    cout << "Average Parallel Edges Insertion time: " << fixed << setprecision(2) << time_cumm / num_insertion_batches * 1000 << " milliseconds" << endl;

    cout << "---------------------------------------END-----------------------------------------" << endl << endl << endl;
    
    return 0;
}
