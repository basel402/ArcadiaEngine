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
    // TODO: Define your skip list node structure and necessary variables
    // Hint: You'll need nodes with multiple forward pointers

public:
    ConcreteLeaderboard() {
        // TODO: Initialize your skip list
    }

    void addScore(int playerID, int score) override {
        // TODO: Implement skip list insertion
        // Remember to maintain descending order by score
    }

    void removePlayer(int playerID) override {
        // TODO: Implement skip list deletion
    }

    vector<int> getTopN(int n) override {
        // TODO: Return top N player IDs in descending score order
        return {};
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


public:
    ConcreteAuctionTree() {
        NIL = new Node();
        NIL->color = BLACK;
        NIL->left = NIL;
        NIL->right = NIL;
        NIL->parent = NIL;
        root = NIL;
    }
    void insertItem(int itemID, int price) override {
        // TODO: Implement Red-Black Tree insertion
        // Remember to maintain RB-Tree properties with rotations and recoloring
    }

    void deleteItem(int itemID) override {
        // TODO: Implement Red-Black Tree deletion
        // This is complex - handle all cases carefully
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
    return 0;
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
    return 0;
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
// =========================================================

bool WorldNavigator::pathExists(int n, vector<vector<int>>& edges, int source, int dest) {
    // TODO: Implement path existence check using BFS or DFS
    // edges are bidirectional
    return false;
}

long long WorldNavigator::minBribeCost(int n, int m, long long goldRate, long long silverRate,
                                       vector<vector<int>>& roadData) {
    // TODO: Implement Minimum Spanning Tree (Kruskal's or Prim's)
    // roadData[i] = {u, v, goldCost, silverCost}
    // Total cost = goldCost * goldRate + silverCost * silverRate
    // Return -1 if graph cannot be fully connected
    return -1;
}

string WorldNavigator::sumMinDistancesBinary(int n, vector<vector<int>>& roads) {
    // TODO: Implement All-Pairs Shortest Path (Floyd-Warshall)
    // Sum all shortest distances between unique pairs (i < j)
    // Return the sum as a binary string
    // Hint: Handle large numbers carefully
    return "0";
}

// =========================================================
// PART D: SERVER KERNEL (Greedy)
// =========================================================

int ServerKernel::minIntervals(vector<char>& tasks, int n) {
    // TODO: Implement task scheduler with cooling time
    // Same task must wait 'n' intervals before running again
    // Return minimum total intervals needed (including idle time)
    // Hint: Use greedy approach with frequency counting
    return 0;
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
