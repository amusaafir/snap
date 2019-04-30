/*TODO: get rid of duplication (undirected/directed)
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Snap.h"
#include <unordered_map>
#include <unordered_set>
#include "stdafx.h"
#include <sys/time.h>

#include <time.h>
int AMOUNT_OF_RANDOM_NODES_CLUSTERING_COEFFICIENT;
int AMOUNT_OF_RANDOM_NODES_DIAMETER;
bool COUNT_NODES_AND_EDGES = true;
bool CALCULATE_AVERAGE_DEGREE = true;
bool CALCULATE_DENSITY = true;
bool CALCULATE_CON_COMP = true;
bool CALCULATE_CLUSTERING_COEFFICIENT = true;
bool CALCULATE_DIAMETER_AND_AVG_PATH = true;

///
//
int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

void timeval_print(struct timeval *tv)
{
    char buffer[30];
    time_t curtime;

    printf("%ld.%06ld", tv->tv_sec, tv->tv_usec);
    curtime = tv->tv_sec;
    strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
    printf(" = %s.%06ld\n", buffer, tv->tv_usec);
}
//
//


int add_vertex_as_coordinate(std::unordered_map<std::string, int>& map_from_edge_to_coordinate, std::string vertex, int* coordinate) {
    if (map_from_edge_to_coordinate.count(vertex)) {
        return map_from_edge_to_coordinate[vertex];
    }
    else {
        map_from_edge_to_coordinate[vertex] = *coordinate;
        return (*coordinate)++;
    }
}

void read_undirected_graph(std::string file_path) {

    std::string inputGraphName = file_path.substr(file_path.find_last_of("/\\") + 1);

    std::cout << inputGraphName<< std::endl;
    const char *array = inputGraphName.c_str();
    TStr inputGraphNameSnapFormat(array);
    std::ofstream outputFile("summary_" +  inputGraphName + ".csv");

    outputFile << "Input graph\t" + inputGraphName + "\n";
    outputFile << "Type\tUndirected\n";


    std::cout << "Undirected graph\t" << std::endl;
    PUNGraph Graph = TUNGraph::New();

    std::unordered_map<std::string, int> map_from_edge_to_coordinate;

    std::ifstream infile(file_path.c_str());
    std::cout << file_path << std::endl;
    std::string vertex_source;
    std::string vertex_destination;
    // std::string edge_weight; // note - this is new
    int current_coordinate = 0;

    std::unordered_set<std::string> vertices;

    int ownEdgeCounter = 0;

    while (infile >> vertex_source >> vertex_destination) {
        //std::cout << "\n" << "-" << vertex_source << "-" << vertex_destination << "-" << std::endl;
        vertices.insert(vertex_source);
        vertices.insert(vertex_destination);
        int iVertex_source = add_vertex_as_coordinate(map_from_edge_to_coordinate, vertex_source, &current_coordinate);
        int iVertex_destination = add_vertex_as_coordinate(map_from_edge_to_coordinate, vertex_destination, &current_coordinate);

        if(iVertex_source<0) {
            printf("Smaller than zero found (source)\t%d", iVertex_source);
        }

        if(iVertex_destination<0) {
            printf("Smaller than zero found (destination)\t%d", iVertex_destination);
        }

        ownEdgeCounter++;

        Graph->AddEdge2(iVertex_source, iVertex_destination);
        //std::cout << "-" << iVertex_source << "," << iVertex_destination <<"-" << std::endl;
    }

    infile.close();

    std::cout << "Done reading the graph." << std::endl;
    std::cout << "Size vertices set\t" << vertices.size() << std::endl;
    std::cout << "Own edge count\t" << ownEdgeCounter << std::endl;

    outputFile << "Vertices\t" + std::to_string(vertices.size()) + "\n";
    outputFile << "Edges\t" + std::to_string(ownEdgeCounter) + "\n";

    // Count nodes and edges
    int amountOfNonZeroNodes;
    int amountEdges;
    if (COUNT_NODES_AND_EDGES) {
        amountOfNonZeroNodes = TSnap::CntNonZNodes(Graph);
        amountEdges = TSnap::CntUniqUndirEdges(Graph);

        std::cout << "Unique amount of (non-zero degree) nodes\t" << amountOfNonZeroNodes << std::endl;
        std::cout << "Unique undirected edges\t" << amountEdges << std::endl;

        outputFile << "Unique vertices\t" + std::to_string(amountOfNonZeroNodes) + "\n";
        outputFile << "Unique edges\t" + std::to_string(amountEdges) + "\n";
    }

    // Average degree
    if (CALCULATE_AVERAGE_DEGREE) {
        if (!COUNT_NODES_AND_EDGES) {
            printf("\nCount nodes/edges is required before calculating average degree.");
            exit(1);
        }
        float avgDegree = 2 * float(amountEdges) / float(amountOfNonZeroNodes);
        std::cout << "Avg deg\t" << avgDegree << std::endl;

        outputFile << "Avg deg\t" + std::to_string(avgDegree) + "\n";
    }

    // Density
    if (CALCULATE_DENSITY) {
        if (!COUNT_NODES_AND_EDGES) {
            printf("\nCount nodes/edges is required before calculating density.");
            exit(1);
        }
        float density = (2 * float(amountEdges)) / (float(amountOfNonZeroNodes) * (float(amountOfNonZeroNodes) - 1));
        std::cout << "Density\t" << density << std::endl;
        outputFile << "Density\t" + std::to_string(density) + "\n";
    }

    TSnap::PlotClustCf(Graph,  inputGraphNameSnapFormat, inputGraphNameSnapFormat);

    TSnap::PlotInDegDistr(Graph, inputGraphNameSnapFormat,  inputGraphNameSnapFormat);

    TSnap::PlotOutDegDistr(Graph, inputGraphNameSnapFormat, inputGraphNameSnapFormat);

    // Connected components
    if (CALCULATE_CON_COMP) {
        TCnComV cc;
        TSnap::GetSccs(Graph, cc);
        std::cout << "(Strongly) Connected components\t" << cc.Len() << std::endl;
        outputFile << "(Strongly) Connected components\t" + std::to_string(cc.Len()) + "\n";

        TCnComV wc;
        TSnap::GetWccs(Graph, wc);
        std::cout << "(Weakly) Connected components\t" << wc.Len() << std::endl;
        outputFile << "(Weakly) Connected components\t" + std::to_string(wc.Len()) + "\n";

    }

    // Clustering coefficient
    if (CALCULATE_CLUSTERING_COEFFICIENT) {
        double clusteringCoefficient =TSnap::GetClustCf(Graph, AMOUNT_OF_RANDOM_NODES_CLUSTERING_COEFFICIENT);
        std::cout << "Clustering coefficient\t" << clusteringCoefficient << std::endl;
        outputFile << "Clustering coefficient\t" + std::to_string(clusteringCoefficient) + "\n";
    }

    // Diameter, avg shortest path
    if (CALCULATE_DIAMETER_AND_AVG_PATH) {

        if (AMOUNT_OF_RANDOM_NODES_DIAMETER == -1) {
            AMOUNT_OF_RANDOM_NODES_DIAMETER = amountOfNonZeroNodes;
        }

        double avgDiam;
        int fullDiam;

        TSnap::PlotShortPathDistr(Graph, inputGraphNameSnapFormat, inputGraphNameSnapFormat, AMOUNT_OF_RANDOM_NODES_DIAMETER, avgDiam, fullDiam);


        outputFile << "Diameter\t" + std::to_string(fullDiam) + "\n";
        outputFile << "Avg path\t" + std::to_string(avgDiam) + "\n";

        //std::cout<<"See diameter plt for diameter data" << std::endl;
        //TSnap::GetBfsEffDiam(Graph, AMOUNT_OF_RANDOM_NODES_DIAMETER, false, effDiameter, fullDiam, avgShortestPathLength);
        //std::cout << "Effdiameter: " << effDiameter << ", " << "Full diameter: " << fullDiam << "," << "Avg shortest path: " << avgShortestPathLength << std::endl;
    }

    // Close summary file
    if (outputFile.is_open()) {
        outputFile.close();
    }
}

void read_directed_graph(std::string file_path) {
    std::string inputGraphName = file_path.substr(file_path.find_last_of("/\\") + 1);

    std::cout << inputGraphName<< std::endl;
    const char *array = inputGraphName.c_str();
    TStr inputGraphNameSnapFormat(array);

    std::ofstream outputFile("summary_" +  inputGraphName + ".csv");

    outputFile << "Input graph\t" + inputGraphName + "\n";
    outputFile << "Type\tDirected\n";

    PNGraph Graph = TNGraph::New();

    std::unordered_map<std::string, int> map_from_edge_to_coordinate;
    std::ifstream infile(file_path.c_str());
    std::cout << file_path << std::endl;
    std::string vertex_source;
    std::string vertex_destination;
    int current_coordinate = 0;

    std::unordered_set<std::string> vertices;

    while (infile >> vertex_source >> vertex_destination) {
        //std::cout << "\n" << vertex_source << "," << vertex_destination << std::endl;
        vertices.insert(vertex_source);
        vertices.insert(vertex_destination);

        int iVertex_source = add_vertex_as_coordinate(map_from_edge_to_coordinate, vertex_source, &current_coordinate);
        int iVertex_destination = add_vertex_as_coordinate(map_from_edge_to_coordinate, vertex_destination, &current_coordinate);

        Graph->AddEdge2(iVertex_source, iVertex_destination);
        //std::cout << iVertex_source << "," << iVertex_destination << std::endl;
    }

    infile.close();

    std::cout << "Done reading the graph." << std::endl;
    std::cout << "Size vertices set\t" << vertices.size() << std::endl;

    outputFile << "Vertices\t" + std::to_string(vertices.size()) + "\n";

    // Count nodes and edges
    int amountOfNonZeroNodes;
    int amountEdges;
    if (COUNT_NODES_AND_EDGES) {
        amountOfNonZeroNodes = TSnap::CntNonZNodes(Graph);
        amountEdges = TSnap::CntUniqDirEdges(Graph);

        std::cout << "Unique amount of (non-zero degree) nodes\t" << amountOfNonZeroNodes << std::endl;
        std::cout << "Unique directed edges\t" << amountEdges << std::endl;
    }

    outputFile << "Unique vertices\t" + std::to_string(amountOfNonZeroNodes) + "\n";

    outputFile << "Unique edges\t" + std::to_string(amountEdges) + "\n";

    // Average degree
    if (CALCULATE_AVERAGE_DEGREE) {
        if (!COUNT_NODES_AND_EDGES) {
            printf("\nCount nodes/edges is required before calculating average degree.");
            exit(1);
        }
        float avgDegree = 2 * float(amountEdges) / float(amountOfNonZeroNodes);
        std::cout << "Avg deg\t" << avgDegree << std::endl;
        outputFile << "Avg deg\t" + std::to_string(avgDegree) + "\n";
    }

    // Density
    if (CALCULATE_DENSITY) {
        if (!COUNT_NODES_AND_EDGES) {
            printf("\nCount nodes/edges is required before calculating density.");
            exit(1);
        }
        float density = float(amountEdges) / (float(amountOfNonZeroNodes) * (float(amountOfNonZeroNodes) - 1));
        std::cout << "Density\t" << density << std::endl;

        outputFile << "Density\t" + std::to_string(density) + "\n";
    }

    TSnap::PlotClustCf(Graph, inputGraphNameSnapFormat, inputGraphNameSnapFormat);

    TSnap::PlotInDegDistr(Graph, inputGraphNameSnapFormat,  inputGraphNameSnapFormat);

    TSnap::PlotOutDegDistr(Graph, inputGraphNameSnapFormat, inputGraphNameSnapFormat);

    // Connected components
    if (CALCULATE_CON_COMP) {
        TCnComV cc;
        TSnap::GetSccs(Graph, cc);
        std::cout << "(Strongly) Connected components\t" << cc.Len() << std::endl;
        outputFile << "(Strongly) Connected components\t" + std::to_string(cc.Len()) + "\n";

        TCnComV wc;
        TSnap::GetWccs(Graph, wc);
        std::cout << "(Weakly) Connected components\t" << wc.Len() << std::endl;
        outputFile << "(Weakly) Connected components\t" + std::to_string(wc.Len()) + "\n";
    }

    // Clustering coefficient
    if (CALCULATE_CLUSTERING_COEFFICIENT) {
        double clusteringCoefficient = TSnap::GetClustCf(Graph, AMOUNT_OF_RANDOM_NODES_CLUSTERING_COEFFICIENT);
        std::cout << "Clustering coefficient\t" << clusteringCoefficient << std::endl;
        outputFile << "Clustering coefficient\t" + std::to_string(clusteringCoefficient) + "\n";
    }

    // Diameter, avg shortest path
    if (CALCULATE_DIAMETER_AND_AVG_PATH) {
        int fullDiam;
        double avgDiam;
        if (AMOUNT_OF_RANDOM_NODES_DIAMETER == -1) {
            AMOUNT_OF_RANDOM_NODES_DIAMETER = amountOfNonZeroNodes;
        }
        TSnap::PlotShortPathDistr(Graph, inputGraphNameSnapFormat, inputGraphNameSnapFormat, AMOUNT_OF_RANDOM_NODES_DIAMETER, avgDiam, fullDiam);

        outputFile << "Diameter\t" + std::to_string(fullDiam) + "\n";
        outputFile << "Avg path\t" + std::to_string(avgDiam) + "\n";

        //std::cout<<"See diameter plt for diameter data" << std::endl;
        //TSnap::GetBfsEffDiam(Graph, AMOUNT_OF_RANDOM_NODES_DIAMETER, true, effDiameter, fullDiam, avgShortestPathLength);
        //std::cout << "Effdiameter: " << effDiameter << ", " << "Full diameter: " << fullDiam << "," << "Avg shortest path: " << avgShortestPathLength << std::endl;
    }

    // Close summary file
    if (outputFile.is_open()) {
        outputFile.close();
    }
}

void set_parameters_false() {
    COUNT_NODES_AND_EDGES = false;
    CALCULATE_AVERAGE_DEGREE = false;
    CALCULATE_DENSITY = false;
    CALCULATE_CON_COMP = false;
    CALCULATE_CLUSTERING_COEFFICIENT = false;
    CALCULATE_DIAMETER_AND_AVG_PATH = false;
}

void get_diameter_parameter(int arg_index, char* argv[]) {
    // Determines whether to calculate the exact diameter and average path length
    sscanf(argv[arg_index], "%d", &AMOUNT_OF_RANDOM_NODES_DIAMETER);
    if (strcmp(argv[arg_index], "-1") == 0) {
        std::cout << "(Compute full diameter and average path length enabled)" << std::endl;
        AMOUNT_OF_RANDOM_NODES_DIAMETER = -1;
    }
    else {
        std::cout << "- Compute diameter and average path length from " << AMOUNT_OF_RANDOM_NODES_DIAMETER << " random nodes" << std::endl;
    }
}

void get_cc_parameter(int arg_index, char* argv[]) {
    // Determines whether to calculate the exact clustering coefficient
    sscanf(argv[arg_index], "%d", &AMOUNT_OF_RANDOM_NODES_CLUSTERING_COEFFICIENT);
    if (strcmp(argv[arg_index], "-1") == 0) {
        std::cout << " - Compute exact clustering coefficient enabled" << std::endl;
        AMOUNT_OF_RANDOM_NODES_CLUSTERING_COEFFICIENT = -1;
    }
    else {
        std::cout << "- Compute exact clustering coefficient " << AMOUNT_OF_RANDOM_NODES_CLUSTERING_COEFFICIENT << " random nodes" << std::endl;
    }
}

void start_calculating(char* argv[]) {
    // Start computing for directed/undirected
    if (strcmp(argv[2], "directed") == 0) {
        std::cout << "Computing statistics for DIRECTED graph:" << std::endl;
        read_directed_graph(argv[1]);
    }
    else {
        std::cout << "Computing statistics for UNDIRECTED graph:" << std::endl;
        read_undirected_graph(argv[1]);
    }
}

int main(int argc, char* argv[]) {




    // Determines whether to calculate all or a single property
    if (strcmp(argv[3], "avgdegree") == 0) {
        set_parameters_false();
        COUNT_NODES_AND_EDGES = true;
        CALCULATE_AVERAGE_DEGREE = true;
    }
    else if (strcmp(argv[3], "density") == 0) {
        set_parameters_false();
        COUNT_NODES_AND_EDGES = true;
        CALCULATE_DENSITY = true;
    }
    else if (strcmp(argv[3], "cc") == 0) {
        set_parameters_false();
        CALCULATE_CLUSTERING_COEFFICIENT = true;
        get_cc_parameter(4, argv);
    }
    else if (strcmp(argv[3], "concom") == 0) {
        set_parameters_false();
        CALCULATE_CON_COMP = true;
    }
    else if (strcmp(argv[3], "diameter") == 0) {
        set_parameters_false();
        CALCULATE_DIAMETER_AND_AVG_PATH = true;
        get_diameter_parameter(4, argv);
    }
    else if (strcmp(argv[3], "all") == 0) {
        get_cc_parameter(4, argv);
        get_diameter_parameter(5, argv);
    }
    else {
        printf("\nIncorrect parameter specification.");
        exit(1);
    }

    start_calculating(argv);

    std::cout << "\nFinished" << std::endl;


    return 0;
}

