

/**
 * @file tvb.cpp
 * @brief Implementation of TVB (Tuple-Valued Bitmask) algorithm for dynamic graph processing
 *
 * This code implements the TVB data structure with its two algorithms for maintaining a MIS
 * in a dynamic graph. It includes parallel implementations for both insertion (PDMI_TVB)
 * and deletion (PDMD_TVB) algorithms.
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

// Number of threads to use for parallel processing
int n_threads = 75;

using namespace std;

// Vector to store execution time for each thread
vector<double> time_per_thread;

/**
 * @class TVB
 * @brief Represents a Tuple-Value Bitmask (TVB) graph structure
 * 
 * This class implements the TVB graph data structure 
 */
class TVB
{
public:
    std::vector<int> edges;   // stores destination nodes of edges
    std::vector<int> indices; // stores index in 'edges' for the start of each node's edge list
    std::vector<bool> membership;
    int maxNode = -1;


    /**
     * @brief Checks if a given set of nodes forms an independent set
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
     * @return true if nodes are adjacent, false otherwise
     */
    bool isAdjacent(int u, int v) const
    {
        if (u >= indices.size() - 1 || v >= indices.size() - 1)
            return false;

        for (int i = indices[u]; i < indices[u + 1]; ++i)
        {
            if (edges[i] == v)
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
        if (u < indices.size() - 1)
        {
            for (int i = indices[u]; i < indices[u + 1]; ++i)
            {
                neighbors.push_back(edges[i]);
            }
        }
        if (u == indices.size() - 1)
        {
            for (int i = indices[u]; i < edges.size(); i++)
            {
                neighbors.push_back(edges[i]);
            }
        }
        return neighbors;
    }
};

/**
 * @brief Reads a Compressed Sparse Row (CSR) representation from a file
 * @param filename Name of the file to read from
 * @return TVB object with the graph data
 */
TVB readCSR(string filename)
{
    TVB tvb;
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
                    edges.push_back(std::stoi(line)+1);
                }
                break;
            }
        }
    }
    // Initialize TVB object with read data
    tvb.maxNode = indices.size();
    tvb.indices.resize(indices.size() + 1, 0);
    tvb.edges.resize(edges.size(), 0);
    tvb.membership.resize(tvb.maxNode + 1);

    tvb.indices[0] = 0;
    for (int i = 1; i <= indices.size(); i++)
    {
        tvb.indices[i] = indices[i - 1];
    }
    for (int i = 0; i < edges.size(); i++)
    {
        tvb.edges[i] = edges[i] ;
    }
    cout << "Node->" << tvb.maxNode + 1 << endl;
    cout << "Edges->" << tvb.edges.size() << endl;

    return tvb;
}

/**
 * @brief Performs Parallel Dynamic MIS Insertion (PDMI) algorithm on TVB
 * @param tvb TVB object representing the graph
 * @param edgeListBatch Batch of edges to process
 */
void PDMI_TVB(TVB &tvb,vector<pair<int, int>> &edgeListBatch)
{
    int block, tid, start, end, it;
    omp_set_num_threads(n_threads);
    double startTime, endTime;
    block = ceil(static_cast<float>(edgeListBatch.size()) / n_threads);
#pragma omp parallel private(tid, start, end, it) shared(tvb)
    {
        startTime = omp_get_wtime();
        tid = omp_get_thread_num();
        start = tid * block;
        end = std::min((tid + 1) * block, static_cast<int>(edgeListBatch.size()));
        for (int j = start; j < end; j++)
        {
            if (tvb.membership[edgeListBatch[j].first] == true && tvb.membership[edgeListBatch[j].second] == true)
            {
                tvb.membership[edgeListBatch[j].first] = false;
                vector<int> neighbours=tvb.getNeighbors(edgeListBatch[j].first);
                for(int k=0;k<neighbours.size();k++)
                {
                bool shouldWeAdd=true;
                    vector<int> nofn=tvb.getNeighbors(neighbours[k]);
                    for(int l=0;l<tvb.getNeighbors(neighbours[k]).size();l++){
                        if(tvb.membership[nofn[l]]){
                            shouldWeAdd=false;
                            break;
                        }
                    }
                if(shouldWeAdd){
                    tvb.membership[neighbours[k]]=true;
                }
                }
            }
        }
        endTime = omp_get_wtime();
        time_per_thread[tid] = (endTime - startTime);
    }
}

/**
 * @brief Performs Parallel Dynamic MIS Deletion (PDMD) algorithm on TVB
 * @param tvb TVB object representing the graph
 * @param edgeListBatch Batch of edges to process
 */
void PDMD_TVB(TVB &tvb, vector<pair<int, int>> &edgeListBatch)
{
    int block, tid, start, end, it;
    double startTime = 0, endTime = 0;
    block = ceil(static_cast<float>(edgeListBatch.size()) / n_threads);
    omp_set_num_threads(n_threads);
#pragma omp parallel private(tid, start, end, it) shared(tvb)
    {
        startTime = omp_get_wtime();
        tid = omp_get_thread_num();
        start = tid * block;
        end = std::min((tid + 1) * block, static_cast<int>(edgeListBatch.size()));

        for (int i = start; i < end; i++)
        {
            pair<int, int> edge = edgeListBatch[i];
            bool firstIn = false, secondIn = false;
            firstIn = tvb.membership[edge.first];
            secondIn = tvb.membership[edge.second];
            if (firstIn == true && secondIn == true)
            {
                continue;
            }
            if (firstIn == true)
            {
                vector<int> neighbours = tvb.getNeighbors(edge.second);
                bool shouldWeAdd = true;

                for (int i = 0; i < neighbours.size(); i++)
                {
                    if (tvb.membership[neighbours[i]] == true && neighbours[i] != edge.first)
                    {
                        shouldWeAdd = false;
                    }
                }

                if (shouldWeAdd && tvb.membership[edge.second] == false)
                {
                    tvb.membership[edge.second] = true;
                }
            }
            else if (secondIn == true)
            {
                vector<int> neighbours = tvb.getNeighbors(edge.first);
                bool shouldWeAdd = true;

                for (int i = 0; i < neighbours.size(); i++)
                {
                    if (tvb.membership[neighbours[i]] == true && neighbours[i] != edge.second)
                    {
                        shouldWeAdd = false;
                    }
                }

                if (shouldWeAdd && tvb.membership[edge.first] == false)
                {
                    tvb.membership[edge.first] = true;

                }
            }
        }
        endTime = omp_get_wtime();
        time_per_thread[tid] = endTime - startTime;
    }
}

/**
 * @brief Main function to run the TVB algorithm
 * @param argc Number of command-line arguments
 * @param argv Array of command-line arguments
 * @return 0 
 */

int main(int argc, char *argv[])
{
    if (argc < 7)
    {
        std::cerr << "Usage: " << argv[0] << " <CSR> <MIS> <Deletion_folder> <Insertion_folder> <num_threads> <Num_batch for deletion> <Num_batch for insertion>" << std::endl;
        return 1;
    }
    cout<<"----------------------------------------Start-----------------------------------------"<<endl;
    cout<<"----------------------------------------TVB-----------------------------------------"<<endl;

    std::string filename1 = argv[1];
    TVB tvb;
    tvb = readCSR(filename1);
    std::vector<int> maximalIndependentSet;
    std::string line;
    // Read MIS---------------------------------------------------------
    std::ifstream inputFile(argv[2]);
    if (inputFile.is_open())
    {
        while (std::getline(inputFile, line))
        {
            maximalIndependentSet.push_back(std::stoi(line)+1);
            tvb.membership[std::stoi(line)+1] = true;
        }
        inputFile.close();
    }
    else
    {
        std::cerr << "Unable to open the file for reading" << std::endl;
    }
    // Reading Batches---------------------------------------------------------
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


 // Print experiment setup information
    int initialCard = maximalIndependentSet.size();
    n_threads = stoi(argv[5]);
    double startTime = 0, endTime = 0, time_cumm = 0;
    int num_deletion_batches = std::stoi(argv[6]);
    int num_insertion_batches = std::stoi(argv[7]);
    cout<< "Cores USed-: "<< n_threads<<endl;
    cout << "TVB Name->" << filename1 << endl;
    cout <<"Batch Size->"<<edgeListsForInsertion[0].size() <<" | Number of batches->"<<num_insertion_batches<< " | Incremental Edges->" << edgeListsForInsertion[0].size() * num_insertion_batches << endl;
    cout <<"Batch Size->"<<edgeListsForDeletion[0].size() <<" | Number of batches->"<<num_deletion_batches<< " | Deletion Edges->" << edgeListsForInsertion[0].size() * num_insertion_batches << endl;
    cout << "MIS set size->" << maximalIndependentSet.size() << endl;


//----------------------------------------------Insertion---------------------------------------------------------

    auto start = omp_get_wtime();
    double time_cum = 0;
    time_per_thread.resize(n_threads, 0.0);
    int counter=0;
    for (int i = 0; i < num_insertion_batches; i++)
    {
        PDMI_TVB(tvb,edgeListsForInsertion[i]);
        double maxi = 0;
        for (int k = 0; k < time_per_thread.size(); k++)
        {
            // calculate the max element
            if (maxi < time_per_thread[k])
            {
                maxi = time_per_thread[k];
            }
        }
        time_cum += maxi;
    }
    auto end = omp_get_wtime();
    int c = 0,p=0;
    for (int i = 0; i < tvb.membership.size(); i++)
    {
        if (tvb.membership[i])
        {
            c++;
        }
    }

    cout << endl;
    cout<< "Initial Card"<<initialCard<<endl;
    cout << "-------------------------------PDMI-TVB------------------------------" << endl;
    cout << "Cardinality changed by->(Insertion)" << abs(initialCard - c) << endl;
    cout << "Total Incremental Edges Insertion time: " << fixed << setprecision(5) << time_cum * 1000 << " milliseconds" << endl;
    cout << "Average Incremental Edges Insertion time: " << fixed << setprecision(5) << time_cum / num_insertion_batches * 1000 << " milliseconds" << endl;
    cout<<endl;
    time_per_thread.clear();
    time_per_thread.resize(n_threads, 0.0);
    time_cum = 0;

    //----------------------------------------------Deletion---------------------------------------------------------

    start = omp_get_wtime();
    for (int i = 0; i < num_deletion_batches; i++)
    {
        PDMD_TVB(tvb, edgeListsForDeletion[i]);
        double maxi = 0;
        for (int k = 0; k < time_per_thread.size(); k++)
        {
            // calculate the max element
            if (maxi < time_per_thread[k])
            {
                maxi = time_per_thread[k];
            }
        }
        time_cum += maxi;
    }
    end = omp_get_wtime();
    int d = 0;
        for (int i = 0; i < tvb.membership.size(); i++)
    {
        if (tvb.membership[i])
        {
            d++;
        }
    }

    cout << "-------------------------------PDMD-TVBL------------------------------" << endl;
    cout << "Cardinality changed by->(Deletion)" << abs(d -c ) << endl;
    cout << "Total Incremental Edges Deletion time: " << fixed << setprecision(5) << time_cum * 1000 << " milliseconds" << endl;
    cout << "Average Incremental Edges Deletion time: " << fixed << setprecision(5) << time_cum / num_deletion_batches * 1000 << " milliseconds" << endl;
    cout << "---------------------------------------END-----------------------------------------" << endl;

}
