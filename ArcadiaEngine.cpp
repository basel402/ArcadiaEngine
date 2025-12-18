// ArcadiaEngine.cpp - STUDENT TEMPLATE
// TODO: Implement all the functions below according to the assignment requirements

#include "ArcadiaEngine.h"
#include <algorithm>
#include <queue>
#include <numeric>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>
#include <cstdint>


using namespace std;

// =========================================================
// PART A: DATA STRUCTURES (Concrete Implementations)
// =========================================================


class ConcretePlayerTable : public PlayerTable {
private:
    #define TABLE_SIZE 101
    vector<int> keys;        // -1 means empty
    vector<string> values;

    int bitsNeeded() const {
        int bits = 0;
        int n = TABLE_SIZE - 1;
        while (n > 0) {
            bits++;
            n >>= 1;
        }
        return bits == 0 ? 1 : bits;
    }

    // Mid-Square hash
    int hash1MidSquare(int key) const {
        uint64_t k = static_cast<uint64_t>(key);
        uint64_t sq = k * k;
        int r = bitsNeeded();
        int shift = (64 - r) / 2;
        uint64_t middle = (sq >> shift) & ((1ULL << r) - 1);
        return static_cast<int>(middle % TABLE_SIZE);
    }

    // Second hash (step size)
    int hash2(int key) const {
        return 1 + ((key >= 0 ? key : -key) % (TABLE_SIZE - 1));
    }

public:
    ConcretePlayerTable()
        : keys(TABLE_SIZE, -1), values(TABLE_SIZE)
    {}

    void insert(int playerID, string name) override {
        if (playerID < 0) return;

        int h1 = hash1MidSquare(playerID);
        int h2 = hash2(playerID);

        for (int i = 0; i < TABLE_SIZE; i++) {
            int idx = (h1 + i * h2) % TABLE_SIZE;

            if (keys[idx] == -1) {
                keys[idx] = playerID;
                values[idx] = name;
                return;
            }

            if (keys[idx] == playerID) {
                values[idx] = name;
                return;
            }
        }

        // table completely full
        cout << "Table is Full" << endl;
    }

    string search(int playerID) override {
        if (playerID < 0) return "";

        int h1 = hash1MidSquare(playerID);
        int h2 = hash2(playerID);

        for (int i = 0; i < TABLE_SIZE; i++) {
            int idx = (h1 + i * h2) % TABLE_SIZE;

            if (keys[idx] == -1) {
                return "";
            }

            if (keys[idx] == playerID) {
                return values[idx];
            }
        }

        return "";
    }
};

// --- 2. Leaderboard (Skip List) ---

class ConcreteLeaderboard : public Leaderboard {
private:
    struct Node {
        int playerID;
        int score;
        vector<Node*> forward;

        Node(int lvl, int id, int sc)
            : playerID(id), score(sc), forward(lvl + 1, nullptr) {}
    };

    const int MAX_LEVEL = 16;
    const double P = 0.5;

    int currentLevel;
    Node* head;

    int randomLevel() {
        int lvl = 0;
        while (((double)rand() / RAND_MAX) < P && lvl < MAX_LEVEL)
            lvl++;
        return lvl;
    }

    bool comesBefore(int id1, int sc1, int id2, int sc2) {
        if (sc1 != sc2)
            return sc1 > sc2;     // descending score
        return id1 < id2;         // ascending ID
    }

public:
    ConcreteLeaderboard() {
        currentLevel = 0;
        head = new Node(MAX_LEVEL, -1, -1);
    }

    void addScore(int playerID, int score) override {
        vector<Node*> update(MAX_LEVEL + 1, nullptr);
        Node* curr = head;

        for (int i = currentLevel; i >= 0; i--) {
            while (curr->forward[i] &&
                   comesBefore(curr->forward[i]->playerID,
                               curr->forward[i]->score,
                               playerID, score)) {
                curr = curr->forward[i];
            }
            update[i] = curr;
        }

        int lvl = randomLevel();
        if (lvl > currentLevel) {
            for (int i = currentLevel + 1; i <= lvl; i++)
                update[i] = head;
            currentLevel = lvl;
        }

        Node* newNode = new Node(lvl, playerID, score);
        for (int i = 0; i <= lvl; i++) {
            newNode->forward[i] = update[i]->forward[i];
            update[i]->forward[i] = newNode;
        }
    }

    void removePlayer(int playerID) override {
        Node* curr = head;
        Node* target = nullptr;

        while (curr->forward[0]) {
            if (curr->forward[0]->playerID == playerID) {
                target = curr->forward[0];
                break;
            }
            curr = curr->forward[0];
        }

        if (!target) return;

        vector<Node*> update(MAX_LEVEL + 1, nullptr);
        curr = head;

        for (int i = currentLevel; i >= 0; i--) {
            while (curr->forward[i] && curr->forward[i] != target)
                curr = curr->forward[i];
            update[i] = curr;
        }

        for (int i = 0; i <= currentLevel; i++) {
            if (update[i]->forward[i] == target)
                update[i]->forward[i] = target->forward[i];
        }

        delete target;

        while (currentLevel > 0 && head->forward[currentLevel] == nullptr)
            currentLevel--;
    }

    vector<int> getTopN(int n) override {
        vector<int> result;
        Node* curr = head->forward[0];

        while (curr && n--) {
            result.push_back(curr->playerID);
            curr = curr->forward[0];
        }

        return result;
    }
};

// --- 3. AuctionTree (Red-Black Tree) ---

class ConcreteAuctionTree : public AuctionTree {
private:
    enum Color { RED, BLACK };

    struct Node {
        int itemID;
        int price;
        Color color;
        Node* left;
        Node* right;
        Node* parent;
    };

    Node* root;
    Node* NIL;

    // ------------------ Node Creation ------------------
    Node* createNode(int itemID, int price) {
        Node* n = new Node();
        n->itemID = itemID;
        n->price = price;
        n->color = RED;
        n->left = NIL;
        n->right = NIL;
        n->parent = NIL;
        return n;
    }

    // ------------------ Rotations ------------------
    void leftRotate(Node* x) {
        Node* y = x->right;
        x->right = y->left;
        if (y->left != NIL)
            y->left->parent = x;

        y->parent = x->parent;
        if (x->parent == NIL)
            root = y;
        else if (x == x->parent->left)
            x->parent->left = y;
        else
            x->parent->right = y;

        y->left = x;
        x->parent = y;
    }

    void rightRotate(Node* y) {
        Node* x = y->left;
        y->left = x->right;
        if (x->right != NIL)
            x->right->parent = y;

        x->parent = y->parent;
        if (y->parent == NIL)
            root = x;
        else if (y == y->parent->right)
            y->parent->right = x;
        else
            y->parent->left = x;

        x->right = y;
        y->parent = x;
    }

    // ------------------ Helpers for Deletion ------------------
    Node* searchByID(Node* node, int id) {
        if (node == NIL) return nullptr;
        if (node->itemID == id) return node;

        Node* leftResult = searchByID(node->left, id);
        if (leftResult != nullptr) return leftResult;

        return searchByID(node->right, id);
    }

    Node* minimum(Node* node) {
        while (node->left != NIL)
            node = node->left;
        return node;
    }

    void rbTransplant(Node* u, Node* v) {
        if (u->parent == NIL)
            root = v;
        else if (u == u->parent->left)
            u->parent->left = v;
        else
            u->parent->right = v;

        v->parent = u->parent;
    }

    // ------------------ Delete Fixup ------------------
    void rbDeleteFixup(Node* x) {
        while (x != root && x->color == BLACK) {
            if (x == x->parent->left) {
                Node* w = x->parent->right;

                if (w->color == RED) {
                    w->color = BLACK;
                    x->parent->color = RED;
                    leftRotate(x->parent);
                    w = x->parent->right;
                }

                if (w->left->color == BLACK && w->right->color == BLACK) {
                    w->color = RED;
                    x = x->parent;
                } else {
                    if (w->right->color == BLACK) {
                        w->left->color = BLACK;
                        w->color = RED;
                        rightRotate(w);
                        w = x->parent->right;
                    }

                    w->color = x->parent->color;
                    x->parent->color = BLACK;
                    w->right->color = BLACK;
                    leftRotate(x->parent);
                    x = root;
                }
            } else {
                Node* w = x->parent->left;

                if (w->color == RED) {
                    w->color = BLACK;
                    x->parent->color = RED;
                    rightRotate(x->parent);
                    w = x->parent->left;
                }

                if (w->right->color == BLACK && w->left->color == BLACK) {
                    w->color = RED;
                    x = x->parent;
                } else {
                    if (w->left->color == BLACK) {
                        w->right->color = BLACK;
                        w->color = RED;
                        leftRotate(w);
                        w = x->parent->left;
                    }

                    w->color = x->parent->color;
                    x->parent->color = BLACK;
                    w->left->color = BLACK;
                    rightRotate(x->parent);
                    x = root;
                }
            }
        }
        x->color = BLACK;
    }

public:
    // ------------------ Constructor ------------------
    ConcreteAuctionTree() {
        NIL = new Node();
        NIL->color = BLACK;
        NIL->left = NIL;
        NIL->right = NIL;
        NIL->parent = NIL;
        root = NIL;
    }

    // ------------------ Insert (Price, ID) ------------------
    void insertItem(int itemID, int price) override {
        Node* z = createNode(itemID, price);
        Node* y = NIL;
        Node* x = root;

        while (x != NIL) {
            y = x;
            if (z->price < x->price ||
               (z->price == x->price && z->itemID < x->itemID))
                x = x->left;
            else
                x = x->right;
        }

        z->parent = y;
        if (y == NIL)
            root = z;
        else if (z->price < y->price ||
                (z->price == y->price && z->itemID < y->itemID))
            y->left = z;
        else
            y->right = z;

        z->left = NIL;
        z->right = NIL;
        z->color = RED;

        while (z->parent->color == RED) {
            if (z->parent == z->parent->parent->left) {
                Node* u = z->parent->parent->right;
                if (u->color == RED) {
                    z->parent->color = BLACK;
                    u->color = BLACK;
                    z->parent->parent->color = RED;
                    z = z->parent->parent;
                } else {
                    if (z == z->parent->right) {
                        z = z->parent;
                        leftRotate(z);
                    }
                    z->parent->color = BLACK;
                    z->parent->parent->color = RED;
                    rightRotate(z->parent->parent);
                }
            } else {
                Node* u = z->parent->parent->left;
                if (u->color == RED) {
                    z->parent->color = BLACK;
                    u->color = BLACK;
                    z->parent->parent->color = RED;
                    z = z->parent->parent;
                } else {
                    if (z == z->parent->left) {
                        z = z->parent;
                        rightRotate(z);
                    }
                    z->parent->color = BLACK;
                    z->parent->parent->color = RED;
                    leftRotate(z->parent->parent);
                }
            }
        }
        root->color = BLACK;
    }

    // ------------------ Delete by Item ID ------------------
    void deleteItem(int itemID) override {
        Node* z = searchByID(root, itemID);
        if (z == nullptr || z == NIL) return;

        Node* y = z;
        Node* x;
        Color yOriginalColor = y->color;

        if (z->left == NIL) {
            x = z->right;
            rbTransplant(z, z->right);
        }
        else if (z->right == NIL) {
            x = z->left;
            rbTransplant(z, z->left);
        }
        else {
            y = minimum(z->right);
            yOriginalColor = y->color;
            x = y->right;

            if (y->parent == z) {
                x->parent = y;
            } else {
                rbTransplant(y, y->right);
                y->right = z->right;
                y->right->parent = y;
            }

            rbTransplant(z, y);
            y->left = z->left;
            y->left->parent = y;
            y->color = z->color;
        }

        if (yOriginalColor == BLACK)
            rbDeleteFixup(x);

        delete z;
    }
};

// =========================================================
// PART B: INVENTORY SYSTEM (Dynamic Programming)
// =========================================================

int InventorySystem::optimizeLootSplit(int n, vector<int>& coins) {
    long long totalSum = 0;
    for (int coin : coins) {
        totalSum += coin;
    }

    int target = totalSum / 2;

    vector<bool> dp(target + 1, false);
    dp[0] = true;

    for (int coin : coins) {
        for (int j = target; j >= coin; j--) {
            if (dp[j - coin]) {
                dp[j] = true;
            }
        }
    }

    int bestSum = 0;
    for (int i = target; i >= 0; i--) {
        if (dp[i]) {
            bestSum = i;
            break;
        }
    }

    return (int)(totalSum - 2 * bestSum);
}

int InventorySystem::maximizeCarryValue(int capacity, vector<pair<int, int>>& items) {
    vector<int> dp(capacity + 1, 0);

    for (const auto& item : items) {
        int weight = item.first;
        int value = item.second;

        for (int w = capacity; w >= weight; w--) {
            dp[w] = max(dp[w], dp[w - weight] + value);
        }
    }

    return dp[capacity];
}

long long InventorySystem::countStringPossibilities(string s) {
    const long long MOD = 1000000007;
    int n = s.length();

    vector<long long> dp(n + 1, 0);
    dp[0] = 1;
    dp[1] = 1;

    for (int i = 2; i <= n; i++) {
        dp[i] = dp[i - 1];
        if (s[i - 1] == 'u' && s[i - 2] == 'u') dp[i] = (dp[i] + dp[i - 2]) % MOD;
        if (s[i - 1] == 'n' && s[i - 2] == 'n') dp[i] = (dp[i] + dp[i - 2]) % MOD;
    }

    return dp[n]%MOD;
}

// =========================================================
// PART C: WORLD NAVIGATOR (Graphs)
// 1) Safe Passage: Path Existence (UNDIRECTED graph)
bool WorldNavigator::pathExists(int n, vector<vector<int>>& edges, int source, int dest) {
    if (n == 0) return false;
    if (source < 0 || source >= n) return false;
    if (dest < 0 || dest >= n) return false;
    if (source == dest) return true;

    vector<vector<int>> adj(n);
    for (auto &e : edges) {
        adj[e[0]].push_back(e[1]);
        adj[e[1]].push_back(e[0]);
    }

    vector<bool> visited(n, false);
    queue<int> q;
    q.push(source);
    visited[source] = true;

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) {
            if (!visited[v]) {
                if (v == dest) return true;
                visited[v] = true;
                q.push(v);
            }
        }
    }
    return false;
}

// 2) The Bribe: Minimum Spanning Tree (Prim's Algorithm)
long long WorldNavigator::minBribeCost(int n, int m, long long goldRate, long long silverRate,
                                      vector<vector<int>>& roadData) {
    vector<vector<pair<int,long long>>> adj(n);

    for (auto &r : roadData) {
        int u = r[0];
        int v = r[1];
        long long cost = r[2] * goldRate + r[3] * silverRate;
        adj[u].push_back({v, cost});
        adj[v].push_back({u, cost});
    }

    vector<bool> inMST(n, false);
    priority_queue<pair<long long,int>, vector<pair<long long,int>>, greater<>> pq;

    pq.push({0, 0});
    long long totalCost = 0;
    int visitedCount = 0;

    while (!pq.empty() && visitedCount < n) {
        auto [cost, u] = pq.top(); pq.pop();
        if (inMST[u]) continue;

        inMST[u] = true;
        totalCost += cost;
        visitedCount++;

        for (auto &edge : adj[u]) {
            if (!inMST[edge.first]) {
                pq.push({edge.second, edge.first});
            }
        }
    }

    return (visitedCount == n) ? totalCost : -1;
}

// 3) Teleporter Network: All-Pairs Shortest Path (Floyd–Warshall)
string WorldNavigator::sumMinDistancesBinary(int n, vector<vector<int>>& roads) {
    const long long INF = 1e18;
    vector<vector<long long>> dist(n, vector<long long>(n, INF));

    for (int i = 0; i < n; i++) dist[i][i] = 0;

    for (auto &r : roads) {
        int u = r[0];
        int v = r[1];
        long long w = r[2];
        dist[u][v] = min(dist[u][v], w);
        dist[v][u] = min(dist[v][u], w);
    }

    for (int k = 0; k < n; k++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (dist[i][k] < INF && dist[k][j] < INF)
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);

    long long sum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (dist[i][j] < INF)
                sum += dist[i][j];
        }
    }

    if (sum == 0) return "0";

    string binary = "";
    while (sum > 0) {
        binary = char('0' + (sum % 2)) + binary;
        sum /= 2;
    }
    return binary;
}

// ================= PART D: SERVER KERNEL =================

int ServerKernel::minIntervals(vector<char>& tasks, int n) {
    if (tasks.empty()) return 0;

    vector<int> freq(26, 0);
    for (char c : tasks)
        freq[c - 'A']++;

    int maxFreq = *max_element(freq.begin(), freq.end());
    int countMax = 0;
    for (int f : freq)
        if (f == maxFreq) countMax++;

    int intervals = (maxFreq - 1) * (n + 1) + countMax;
    return max((int)tasks.size(), intervals);
}


// =========================================================
// FACTORY FUNCTIONS (Required for Testing)
// =========================================================

extern "C" {
    PlayerTable* createPlayerTable() { 
        return new ConcretePlayerTable(); 
    }

    Leaderboard* createLeaderboard() { 
        return new ConcreteLeaderboard(); 
    }

    AuctionTree* createAuctionTree() { 
        return new ConcreteAuctionTree(); 
    }
}
