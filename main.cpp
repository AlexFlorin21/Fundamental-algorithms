#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <fstream>
using namespace std;

class Graf {
private:
    unordered_map<int, vector<int>> AdjGraph; // Container pentru reprezentarea listei de adiacenta

    vector<vector<int>> MatrixGraph; // Container pentru stocarea unui graf ca o matrice

    vector<int> Parent; // Vector care tine evidenta parintilor in arbore

    vector<int> Visited;

    // Adauga o muchie intre nodurile x si y in graful dat
    void AddEdge(int x, int y);

    // Realizeaza o sortare topologica a grafului si returneaza ordinea topologica a nodurilor
    vector<int> TopologicalSort();

    void DfsCriticalConnections(int node, int Parent, vector<int>& discovery, vector<int>& minimum, int& time, vector<vector<int>>& criticalConnections);

    // Gaseste si returneaza conexiunile critice din graful dat folosind algoritmul lui Tarjan
    vector<vector<int>> FindCriticalConnections();

    // Verifica daca graful este bipartit
    bool Bipartition();

    bool DfsBipartition(int node, int group, unordered_map<int, int>& groups);

    // Gaseste si returneaza nodurile din graful dat care nu fac parte din cicluri
    bool verifyCycleFreeNodes(int node);

    vector<int> NodesWithoutCycles();

    // Gaseste cea mai scurta distanta de la insula folosind o cautare in latime (BFS)
    int FindShortestDistanceToIsland();

    void dfs(int x, int y, queue<pair<int, int>>& q, vector<vector<int>>& Visited);

    //Gaseste radacina unui nod
    int find(int x);

    //Verifica daca un sistem de ecuatii este valid
    bool IsValidEquations(vector<string>& equations);

public:
    Graf() {
        Parent.resize(26);
        for (int i = 0; i < 26; ++i) {
            Parent[i] = i;
        }
    }

    const unordered_map<int, vector<int>>& getAdjGraph() const {
        return AdjGraph;
    }

    const vector<vector<int>>& getMatrixGraph() const {
        return MatrixGraph;
    }

    void setMatrixGraph(const vector<vector<int>>& InputMatrixGraph) {
        MatrixGraph = InputMatrixGraph;
    }

    void setEdge(int x, int y);
    void menu();
    ~Graf() {}
};
void Graf::setEdge(int x, int y) {

    AddEdge(x, y);

}

void Graf::AddEdge(int x, int y) {

    AdjGraph[x].push_back(y);
    AdjGraph[y].push_back(x);

}


vector<int> Graf::TopologicalSort() {
    vector<int> inDegree(AdjGraph.size(), 0);
    vector<int> courseOrder;
    queue<int> q;

    for (const auto& entry : AdjGraph) {
        for (int course : entry.second) {
            inDegree[course]++;
        }
    }

    for (int course = 0; course < AdjGraph.size(); course++) {
        if (inDegree[course] == 0) {
            q.push(course);
        }
    }

    while (!q.empty()) {
        int course = q.front();
        q.pop();
        courseOrder.push_back(course);

        for (int neighbor : AdjGraph[course]) {
            inDegree[neighbor]--;

            if (inDegree[neighbor] == 0) {
                q.push(neighbor);
            }
        }
    }

    if (courseOrder.size() == AdjGraph.size()) {
        return courseOrder;
    } else {
        return vector<int>();
    }
}


void Graf::DfsCriticalConnections(int node, int Parent, vector<int>& discovery, vector<int>& minimum, int& time, vector<vector<int>>& criticalConnections) {
    discovery[node] = minimum[node] = ++time;

    for (int neighbor : AdjGraph[node]) {
        if (neighbor == Parent) continue;

        if (discovery[neighbor] == 0) {
            DfsCriticalConnections(neighbor, node, discovery, minimum, time, criticalConnections);
            minimum[node] = min(minimum[node], minimum[neighbor]);

            if (minimum[neighbor] > discovery[node]) {
                criticalConnections.push_back({node, neighbor});
            }
        } else {
            minimum[node] = min(minimum[node], discovery[neighbor]);
        }
    }
}

vector<vector<int>> Graf::FindCriticalConnections() {
        vector<vector<int>> criticalConnections;
        int n = AdjGraph.size();
        vector<int> discovery(n, 0);
        vector<int> minimum(n, 0);
        int time = 0;

        for (int i = 0; i < n; i++) {
            if (discovery[i] == 0) {
                DfsCriticalConnections(i, -1, discovery, minimum, time, criticalConnections);
            }
        }

        return criticalConnections;
    }

bool Graf::Bipartition() {
        unordered_map<int, int> groups;

        for (const auto& entry : AdjGraph) {
            int person = entry.first;
            if (groups.find(person) == groups.end() && !DfsBipartition(person, 1, groups)) {
                return false;
            }
        }

        return true;
    }



bool Graf::DfsBipartition(int node, int group, unordered_map<int, int>& groups) {
    if (groups.find(node) != groups.end()) {
        return groups[node] == group;
    }

    groups[node] = group;

    for (int neighborNode : AdjGraph[node]) {
        if (!DfsBipartition(neighborNode, 3 - group, groups)) {
            return false;
        }
    }

    return true;
}

bool Graf:: verifyCycleFreeNodes(int node) {
        if (Visited[node] == 1) {
            return true;
        }

        if (Visited[node] == -1) {
            return false;
        }

        Visited[node] = -1;

        for (int neighbor : MatrixGraph[node]) {
            if (!verifyCycleFreeNodes(neighbor)) {
                return false;
            }
        }

        Visited[node] = 1;
        return true;
    }

vector<int> Graf:: NodesWithoutCycles() {
        int n = MatrixGraph.size();
        vector<int> safeNodes;

        for (int node = 0; node < n; node++) {
            if (verifyCycleFreeNodes(node)) {
                safeNodes.push_back(node);
            }
        }

        sort(safeNodes.begin(), safeNodes.end());

        return safeNodes;
    }


int Graf :: FindShortestDistanceToIsland() {
        int n = MatrixGraph.size();
        vector<vector<int>> visited(n, vector<int>(n, 0));
        queue<pair<int, int>> q;

        bool foundFirstIsland = false;

        for (int i = 0; i < n && !foundFirstIsland; ++i) {
            for (int j = 0; j < n && !foundFirstIsland; ++j) {
                if (MatrixGraph[i][j] == 1) {
                   dfs(i, j, q, visited);
                    foundFirstIsland = true;
                }
            }
        }

        int steps = 0;
        int directions[] = {-1, 0, 1, 0, -1};

        while (!q.empty()) {
            int size = q.size();

            for (int i = 0; i < size; ++i) {
                int x = q.front().first;
                int y = q.front().second;
                q.pop();

                for (int d = 0; d < 4; ++d) {
                    int newX = x + directions[d];
                    int newY = y + directions[d + 1];

                    if (newX >= 0 && newX < n && newY >= 0 && newY < n) {
                        if (visited[newX][newY] == 0) {
                            if (MatrixGraph[newX][newY] == 1) {
                                return steps;
                            }

                            q.push({newX, newY});
                            visited[newX][newY] = 1;
                        }
                    }
                }
            }

            ++steps;
        }

        return -1;
    }


void Graf::dfs(int x, int y, queue<pair<int, int>>& q, vector<vector<int>>& Visited) {
    int n = MatrixGraph.size();

    if (x < 0 || x >= n || y < 0 || y >= n || MatrixGraph[x][y] == 0 || Visited[x][y] == 1) {
        return;
    }

    Visited[x][y] = 1;
    q.push({x, y});

    dfs(x - 1, y, q, Visited);
    dfs(x + 1, y, q, Visited);
    dfs(x, y - 1, q, Visited);
    dfs(x, y + 1, q, Visited);
}


int Graf::find(int x) {
    if (Parent[x] == x) return x;
    return Parent[x] = find(Parent[x]);
}

bool Graf::IsValidEquations(vector<string>& equations) {
    for (const string& equation : equations) {
        if (equation[1] == '=') {
            int x = equation[0] - 'a';
            int y = equation[3] - 'a';
            int rootX = find(x);
            int rootY = find(y);

            if (rootX != rootY) {
                Parent[rootX] = rootY;
            }
        }
    }

    for (const string& equation : equations) {
        if (equation[1] == '!' && find(equation[0] - 'a') == find(equation[3] - 'a')) {
            return false;
        }
    }

    return true;
}


using namespace std;

int main() {
    Graf Graf;

    ifstream inputFile("graph_input.txt");

    int NrNodes, NrEdges;
    inputFile >> NrNodes >> NrEdges;

    for (int i = 0; i < NrEdges; i++) {
        int x, y;
        inputFile >> x >> y;
        Graf.setEdge(x, y);
    }

    inputFile.close();



    return 0;
}
