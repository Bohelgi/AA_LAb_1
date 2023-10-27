#include<iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <locale>

using namespace std;

// Клас для представлення графа
class Graph {
public:
    
    Graph(int numVertices) : numVertices(numVertices) {
        adjacencyMatrix.resize(numVertices, vector<int>(numVertices, 0));
    }

    void addEdge(int v1, int v2, int weight = 1) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = weight;
            adjacencyMatrix[v2][v1] = weight;
        }
    }

    
     void addVertex() {
        numVertices++;
        for (int i = 0; i < numVertices; ++i) {
            adjacencyMatrix[i].push_back(0);
        }
        adjacencyMatrix.push_back(vector<int>(numVertices, 0));
    }

    
    void rmvEdge(int v1, int v2) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = 0;
            adjacencyMatrix[v2][v1] = 0;
        }
    }

    
     void rmvVertex(int v) {
        if (v >= 0 && v < numVertices) {
            numVertices--;
            adjacencyMatrix.erase(adjacencyMatrix.begin() + v);

            for (int i = 0; i < numVertices; ++i) {
                adjacencyMatrix[i].erase(adjacencyMatrix[i].begin() + v);
            }
        }
    }

    // Вивід матриці суміжності графа
    void print() const {
        for (int i = 0; i < numVertices; ++i) {
            for (int j = 0; j < numVertices; ++j) {
                cout << adjacencyMatrix[i][j] << " ";
            }
            cout << endl;
        }
    }

    // Генерація випадкового графа
    void generateRandomGraph() {
        srand(static_cast<unsigned>(time(0)));

        for (int i = 0; i < numVertices; ++i) {
            for (int j = i + 1; j < numVertices; ++j) {
                if (rand() % 2 == 1) {
                    int weight = 1 + rand() % 10;
                    addEdge(i, j, weight);
                }
            }
        }
    }

    // Конвертація у граф у представлення списком списком суміжності
    void convertToListGraph(vector<vector<int>>& adjacencyList) {
        adjacencyList.clear();
        adjacencyList.resize(numVertices);
        for (int i = 0; i < numVertices; ++i) {
            for (int j = 0; j < numVertices; ++j) {
                if (adjacencyMatrix[i][j] != 0) {
                    adjacencyList[i].push_back(j);
                }
            }
        }
    }

protected:
    int numVertices;
    vector<vector<int>> adjacencyMatrix;
};

// Клас для представлення орієнтованго графа
class DirectedGraph : public Graph {
public:
    DirectedGraph(int numVertices) : Graph(numVertices) {}

    void addEdge(int v1, int v2, int weight = 1) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = weight;
        }
    }

    void rmvEdge(int v1, int v2) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = 0;
        }
    }
};

// Клас для представлення неорієнтованго графа
class UndirectedGraph : public Graph {
public:
    UndirectedGraph(int numVertices) : Graph(numVertices) {}

    void rmvVertex(int v) {
        if (v >= 0 && v < numVertices) {
            throw logic_error("Неможливо видалити вершину у неорiєнтованому графу");
        }
    }

    void addVertex() {
        throw logic_error("Неможливо додати вершину до неорiєнтованого графа");
    }
};

// Клас для представлення зваженого графа
class WeightedGraph : public Graph {
public:
    WeightedGraph(int numVertices) : Graph(numVertices) {}

    void addEdge(int v1, int v2, int weight) {
        if (v1 >= 0 && v1 < numVertices && v2 >= 0 && v2 < numVertices) {
            adjacencyMatrix[v1][v2] = weight;
        }
    }


    void generateRandomGraph() {
        Graph::generateRandomGraph();
    }


    // Внутрішні класи для алгоритмів Дейкстри та Беллмана-Форда

    class BellmanFordAlgorithm {
    public:
        BellmanFordAlgorithm(WeightedGraph& graph) : numVertices(graph.getNumVertices()), graph(graph) {}

        // Знаходження найкоротших шляхів від початкової вершини до всіх інших вершин
        void shortestPaths(int source) {
            vector<int> distance(numVertices, INT_MAX);
            distance[source] = 0;

            for (int i = 0; i < numVertices - 1; ++i) {
                for (int v1 = 0; v1 < numVertices; ++v1) {
                    for (int v2 = 0; v2 < numVertices; ++v2) {
                        int weight = graph.getEdgeWeight(v1, v2);
                        if (weight != 0 && distance[v1] != INT_MAX && distance[v1] + weight < distance[v2]) {
                            distance[v2] = distance[v1] + weight;
                        }
                    }
                }
            }
        }

    private:
        int numVertices;
        WeightedGraph& graph;
    };


    
    class DijkstraAlgorithm {
    public:
        DijkstraAlgorithm(WeightedGraph& graph) : numVertices(graph.getNumVertices()), graph(graph) {}

        // Знаходження найкоротших шляхів від початкової вершини до всіх інших вершин
        void shortestPath(int source, int destination) {
            vector<int> distance(numVertices, INT_MAX);
            distance[source] = 0;
            vector<int> visited(numVertices, 0);

            for (int i = 0; i < numVertices - 1; ++i) {
                int minDistance = INT_MAX;
                int minIndex = -1;
                for (int v = 0; v < numVertices; ++v) {
                    if (!visited[v] && distance[v] < minDistance) {
                        minDistance = distance[v];
                        minIndex = v;
                    }
                }

                visited[minIndex] = true;

                for (int v = 0; v < numVertices; ++v) {
                    if (!visited[v] && graph.getEdgeWeight(minIndex, v) != 0 && distance[minIndex] != INT_MAX &&
                        distance[minIndex] + graph.getEdgeWeight(minIndex, v) < distance[v]) {
                        distance[v] = distance[minIndex] + graph.getEdgeWeight(minIndex, v);
                    }
                }
            }
        }

    private:
        int numVertices;
        WeightedGraph& graph;
    };

   
    vector<vector<int>>& getAdjacencyMatrix() {
        return adjacencyMatrix;
    }

    int getNumVertices() const {
        return numVertices;
    }

    int getEdgeWeight(int v1, int v2) const {
        return adjacencyMatrix[v1][v2];
    }
};

int main() {

    setlocale(LC_CTYPE, "ukr");

    cout << "Реалiзацiя графiв" << "\n" <<endl;

    cout << "Орiєнтований граф" << "\n";
    DirectedGraph directedGraph(4);
    directedGraph.generateRandomGraph();
    directedGraph.print();

    cout << "Неорiєнтований граф" << "\n";
    UndirectedGraph undirectedGraph(7);
    undirectedGraph.generateRandomGraph();
    undirectedGraph.print();

    cout << "Зважений граф" << "\n";
    WeightedGraph weightedGraph(5);
    weightedGraph.generateRandomGraph();
    weightedGraph.print();

    cout << "Алгоритм Белмана-Форда" << "\n";
    WeightedGraph::BellmanFordAlgorithm bellmanFordGraph(weightedGraph);
    bellmanFordGraph.shortestPaths(3);

    cout << "Алгоритм Дейкстри" << "\n";
    WeightedGraph::DijkstraAlgorithm dijkstraGraph(weightedGraph);
    dijkstraGraph.shortestPath(0, 3);

    vector<int> graphSizes = { 15, 25, 50, 150, 250, 500 }; // задаю розмірність графам

    for (int size : graphSizes) {
        WeightedGraph weightedGraph(size);
        weightedGraph.generateRandomGraph();

        auto startBellmanFord = chrono::high_resolution_clock::now();
        WeightedGraph::BellmanFordAlgorithm bellmanFordGraph(weightedGraph);
        bellmanFordGraph.shortestPaths(0);
        auto endBellmanFord = chrono::high_resolution_clock::now();
        auto durationBellmanFord = chrono::duration_cast<chrono::microseconds>(endBellmanFord - startBellmanFord);

        auto startDijkstra = chrono::high_resolution_clock::now();
        WeightedGraph::DijkstraAlgorithm dijkstraGraph(weightedGraph);
        dijkstraGraph.shortestPath(0, size - 1);
        auto endDijkstra = chrono::high_resolution_clock::now();
        auto durationDijkstra = chrono::duration_cast<chrono::microseconds>(endDijkstra - startDijkstra);

        cout << "Розмiр графу: " << size << "\n";
        cout << "Час виконання алгоритму Белмана-Форда: " << durationBellmanFord.count() << " мiкросекунди\n";
        cout << "Час виконання алгоритму Дейкстри: " << durationDijkstra.count() << " мiкросекунди\n";
        
    }
        return 0;
}