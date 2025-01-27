/**
 * @file tvbl.cpp
 * @brief Implementation of TVBL (Tuple-Valued Bitmask with Levelling) algorithm for dynamic graph processing
 *
 * This code implements the TVBL data structure with its two algorithms for maintaining a Maximal Independent Set (MIS)
 * in a dynamic graph. It includes parallel implementations for both insertion (PDMI_TVBL)
 * and deletion (PDMD_TVBL) operations.
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

// Number of threads for parallel processing
int n_threads = 75;

using namespace std;

// Vector to store execution time for each thread
vector<double> time_per_thread;

/**
 * @class TVBL
 * @brief Represents the Tuple-Valued Bitmask with Levelling graph structure
 *
 * This class contains the core data structures and methods for the TVBL algorithm,
 * including graph representation, CSR (Compressed Sparse Row) format, and utility functions.
 * Main Variable of focus in this class in the tvbl
 */
class TVBL
{
public:
    std::vector<vector<int>> tvbl;  // Stores TVBL-specific information for each node
    std::vector<int> csredges;      // Stores edges in CSR format 
    std::vector<long int> csrindices;    // Stores indices for CSR format 
    
    int maxNode = -1;  // Stores the maximum node ID in the graph


    /**
     * @brief Checks if a given set of nodes forms an independent set
     * @param mis Vector of node IDs to check
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
     * @brief Checks if two nodes are adjacent
     * @param u First node ID
     * @param v Second node ID
     * @return true if nodes are adjacent, false otherwise
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
     * @param u Node ID
     * @return Vector of neighbor node IDs
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

/**
 * @brief Parallel Dynamic Maximal Independent Set Insertion for TVBL (Insertion)
 * @param tvbl Reference to the TVBL object
 * @param edgeListBatch Vector of edges to be inserted
 */
void PDMI_TVBL(TVBL &tvbl, vector<pair<int, int>> &edgeListBatch)
{
    int block, tid, start, end, it;
    omp_set_num_threads(n_threads);
    double startTime, endTime;
    block = ceil(static_cast<float>(edgeListBatch.size()) / n_threads);

    #pragma omp parallel private(tid, start, end, it) shared(tvbl)
    {
        startTime = omp_get_wtime();
        tid = omp_get_thread_num();
        start = tid * block;
        end = std::min((tid + 1) * block, static_cast<int>(edgeListBatch.size()));

        for (int j = start; j < end; j++)
        {
            if (tvbl.tvbl[edgeListBatch[j].first][1] == true && tvbl.tvbl[edgeListBatch[j].second][1] == true)
            {
                tvbl.tvbl[edgeListBatch[j].first][1] = 0;
                vector<int> neighbours = tvbl.getNeighbors(edgeListBatch[j].first);
                for (int k = 0; k < neighbours.size(); k++)
                {
                    if (tvbl.tvbl[neighbours[k]][2] == 1)
                    {
                        vector<int> nofn = tvbl.getNeighbors(neighbours[k]);
                        tvbl.tvbl[edgeListBatch[j].first][2]++;
                        tvbl.tvbl[neighbours[k]][1] = 1;
                        for(int i = 0; i < nofn.size(); i++)
                        {
                            tvbl.tvbl[nofn[i]][2]++;
                        }
                    }
                    tvbl.tvbl[neighbours[k]][2]--;
                }
            }
        }
        endTime = omp_get_wtime();
        time_per_thread[tid] = (endTime - startTime);
    }
}

/**
 * @brief Parallel Dynamic Maximal Independent Set Deletion for TVBL (Deletion)
 * @param tvbl Reference to the TVBL object
 * @param edgeListBatch Vector of edges to be deleted
 */
void PDMD_TVBL(TVBL &tvbl, vector<pair<int, int>> &edgeListBatch)
{
    int block, tid, start, end, it;
    double startTime = 0, endTime = 0;
    block = ceil(static_cast<float>(edgeListBatch.size()) / n_threads);
    omp_set_num_threads(n_threads);

    #pragma omp parallel private(tid, start, end, it) shared(tvbl)
    {
        startTime = omp_get_wtime();
        tid = omp_get_thread_num();
        start = tid * block;
        end = std::min((tid + 1) * block, static_cast<int>(edgeListBatch.size()));

        for (int l = start; l < end; l++)
        {
            pair<int, int> edge = edgeListBatch[l];
            bool firstIn = false, secondIn = false;
            firstIn = tvbl.tvbl[edge.first][1];
            secondIn = tvbl.tvbl[edge.second][1];

            if (firstIn == true && secondIn == true)
            {
                continue;
            }

            if (firstIn == true)
            {
                if (tvbl.tvbl[edge.second][2] == 1 && tvbl.tvbl[edge.second][1] == false)
                {
                    tvbl.tvbl[edge.second][1] = 1;
                    tvbl.tvbl[edge.second][2] = 0;
                    vector<int> neighbors = tvbl.getNeighbors(edge.second);
                    for(int p = 0; p < neighbors.size(); p++){
                        tvbl.tvbl[neighbors[p]][2]++;
                    }
                }
            }
            else if (secondIn == true)
            {
                if (tvbl.tvbl[edge.first][2] == 1 && tvbl.tvbl[edge.first][1] == false)
                {
                    tvbl.tvbl[edge.first][1] = 1;
                    tvbl.tvbl[edge.first][2] = 0;
                    vector<int> neighbors = tvbl.getNeighbors(edge.first);
                    for(int p = 0; p < neighbors.size(); p++){
                        tvbl.tvbl[neighbors[p]][2]++;
                    }
                }
            }
        }
        endTime = omp_get_wtime();
        time_per_thread[tid] = endTime - startTime;
    }
}

/**
 * @brief Reads CSR format and MIS from files and initializes TVBL object
 * @param filename Path to the CSR file
 * @param misFile Path to the MIS file
 * @return Initialized TVBL object
 */
TVBL readCSR(string filename, string misFile)
{
    TVBL tvbl;
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

    // Initialize TVBL object with read data
    tvbl.maxNode = indices.size();
    tvbl.csrindices.resize(indices.size() + 1, 0);
    tvbl.tvbl.resize(tvbl.maxNode + 1);
    tvbl.csredges.resize(edges.size(), 0);

    for (int i = 0; i < tvbl.maxNode + 1; i++)
    {
        tvbl.tvbl[i].resize(3, 0);
    }

    tvbl.csrindices[0] = 0;
    for (int i = 1; i <= indices.size(); i++)
    {
        tvbl.csrindices[i] = indices[i - 1];
    }

    tvbl.tvbl[0][0] = 0;
    for (int i = 1; i <= indices.size(); i++)
    {
        tvbl.tvbl[i][0] = indices[i - 1];
    }

    for (int i = 0; i < edges.size(); i++)
    {
        tvbl.csredges[i] = edges[i];
    }

    cout << "Nodes->" << tvbl.maxNode + 1 << endl;
    cout << "Edges->" << tvbl.csredges.size() << endl;

    // Read MIS from file
    std::vector<int> maximalIndependentSet;
    std::string line;
    std::ifstream inputFile(misFile);

    if (inputFile.is_open())
    {
        while (std::getline(inputFile, line))
        {
            maximalIndependentSet.push_back(std::stoi(line) + 1);
            tvbl.tvbl[std::stoi(line) + 1][1] = 1;
        }
        inputFile.close();
    }
    else
    {
        std::cerr << "Unable to open the file for reading" << std::endl;
    }

    // Building Levels
    int disj = 0;
    for (int i = 0; i < tvbl.maxNode; i++)
    {
        if (tvbl.tvbl[i][1] == 1)
        {
            tvbl.tvbl[i][2] = 0;
        }
        else
        {
            vector<int> neighbours = tvbl.getNeighbors(i);
            if (neighbours.size() == 0)
            {
                disj++;
            }
            for (int j = 0; j < neighbours.size(); j++)
            {
                if (tvbl.tvbl[neighbours[j]][1] == 1)
                {
                    tvbl.tvbl[i][2]++;
                }
            }
        }
    }

    return tvbl;
}

/**
 * @brief Main function to run TVBL algorithm experiments
 * @param argc Number of command-line arguments
 * @param argv Array of command-line arguments
 * @return 0 on successful execution
 */
int main(int argc, char *argv[])
{
    if (argc < 7)
    {
        std::cerr << "Usage: " << argv[0] << " <CSR> <MIS> <Deletion_folder> <Insertion_folder> <num_threads> <Num_batch for deletion> <Num_batch for insertion>" << std::endl;
        return 1;
    }
    cout << "----------------------------------------Start-----------------------------------------" << endl;
    cout << "----------------------------------------TVBL------------------------------------------" << endl;
    
    std::string filename1 = argv[1];
    std::string filename2 = argv[2];
    TVBL tvbl;

    // Read CSR and MIS files
    tvbl = readCSR(filename1, filename2);
    std::vector<int> maximalIndependentSet;
    std::string line;
    
    // Read MIS
    vector<bool> nodeStats(tvbl.maxNode + 1, false);

    std::ifstream inputFile(argv[2]);
    if (inputFile.is_open())
    {
        while (std::getline(inputFile, line))
        {
            maximalIndependentSet.push_back(std::stoi(line) + 1);
            nodeStats[std::stoi(line) + 1] = true;
        }
        inputFile.close();
    }
    else
    {
        std::cerr << "Unable to open the file for reading" << std::endl;
    }

    // Read Batches for Deletion
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

    // Read Batches for Insertion
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

    int initialCard = maximalIndependentSet.size();
    n_threads = stoi(argv[5]);
    double startTime = 0, endTime = 0, time_cumm = 0;
    int num_deletion_batches = std::stoi(argv[6]);
    int num_insertion_batches = std::stoi(argv[7]);

    // Print experiment setup information
    cout << "Cores Used-: " << n_threads << endl;
    cout << "TVBL Name->" << filename1 << endl;
    cout << "Batch Size->" << edgeListsForInsertion[0].size() << " | Number of batches->" << num_insertion_batches << " | Incremental Edges->" << edgeListsForInsertion[0].size() * num_insertion_batches << endl;
    cout << "Batch Size->" << edgeListsForDeletion[0].size() << " | Number of batches->" << num_deletion_batches << " | Deletion Edges->" << edgeListsForInsertion[0].size() * num_insertion_batches << endl;
    cout << "MIS set size->" << maximalIndependentSet.size() << endl;

    // Insertion Experiment
    auto start = omp_get_wtime();
    double time_cum = 0;
    time_per_thread.resize(n_threads, 0.0);
    int counter = 0;
    TVBL graphforinsertion;
    graphforinsertion = tvbl;
    for (int i = 0; i < num_insertion_batches; i++)
    {
        PDMI_TVBL(graphforinsertion, edgeListsForInsertion[i]);

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

    int tv = 0;
    for (int i = 0; i < graphforinsertion.maxNode; i++)
    {
        if (graphforinsertion.tvbl[i][1] == 1)
        {
            tv++;
        }
    }

    // Print Insertion Results
    cout << endl;
    cout << "Initial Card " << initialCard << endl;
    cout << "-------------------------------PDMI-TVBL------------------------------" << endl;
    cout << "Cardinality changed by->(Insertion): " << abs(initialCard - tv)  << endl;
    cout << "Total Incremental Edges Insertion time: " << fixed << setprecision(5) << time_cum * 1000 << " milliseconds" << endl;
    cout << "Average Incremental Edges Insertion time: " << fixed << setprecision(5) << time_cum / num_insertion_batches * 1000 << " milliseconds" << endl;

    // Deletion Experiment
    time_per_thread.clear();
    time_per_thread.resize(n_threads, 0.0);
    TVBL graphfordeletion;
    graphfordeletion = tvbl;
    time_cum = 0;

    start = omp_get_wtime();
    for (int i = 0; i < num_deletion_batches; i++)
    {
        PDMD_TVBL(graphforinsertion, edgeListsForDeletion[i]);
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

    int c = 0;
    for (int i = 0; i < graphfordeletion.maxNode; i++)
    {
        if (graphforinsertion.tvbl[i][1] == 1)
        {
            c++;
        }
    }

    // Print Deletion Results
    cout << endl;
    cout << "-------------------------------PDMD-TVBL------------------------------" << endl;
    cout << "Cardinality changed by tv->(Deletion): " << abs(c - tv) << endl;
    cout << "Total Incremental Edges Deletion time: " << fixed << setprecision(5) << time_cum * 1000 << " milliseconds" << endl;
    cout << "Average Incremental Edges Deletion time: " << fixed << setprecision(5) << time_cum / num_deletion_batches * 1000 << " milliseconds" << endl;
    
    cout << "---------------------------------------END-----------------------------------------" << endl << endl << endl;

    return 0;
}
