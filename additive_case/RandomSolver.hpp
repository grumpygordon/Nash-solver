#pragma once

#include "Solver.hpp"
#include <bits/stdc++.h>
using namespace std;

// Tries assigning random integers and checks for Nash Equilibrium
struct RandomSolver : Solver {

    pair<int, int> get_score(vector<vector<pair<int, int>>> const &weights, int i1, int j1) {

        vector<int> out(n, -1);
        for (int t = 0; t < 2; t++) {
            int id = vector<int>{i1, j1}[t];
            for (int i = 0; i < (int)assigned[t].size(); i++)
                out[assigned[t][i]] = strategy[t][id][i];
        }
        vector<bool> visited(n);
        int v = start;
        pair<int, int> score = {0, 0};
        while (!visited[v] && v != finish) {
            visited[v] = true;
            int u = out[v];
            assert(u >= 0);
            score.first += weights[v][u].first;
            score.second += weights[v][u].second;
            v = u;
        }
        if (v != finish) {
            pair<int, int> old = score;
            score = {0, 0};
            int fin = v;
            while (true) {
                int u = out[v];
                assert(u >= 0);
                score.first += weights[v][u].first;
                score.second += weights[v][u].second;
                v = u;
                if (v == fin)
                    break;
            }
            if (score.first < 0)
                score.first = numeric_limits<int>::min();
            else if (score.first > 0)
                score.first = numeric_limits<int>::max();
            else
                score.first = old.first;
            if (score.second < 0)
                score.second = numeric_limits<int>::min();
            else if (score.second > 0)
                score.second = numeric_limits<int>::max();
            else
                score.second = old.second;
        }
        return score;
    }

    // Time - number of seconds to run
    void solve(int Time = 600000, int random_seed = -1) {
        if (random_seed == -1)
            random_seed = chrono::high_resolution_clock::now().time_since_epoch().count();
        mt19937 rnd(random_seed);

        clock_t START_TIME = clock();
        int test_count = 0;
        int total_minutes = 0;
        while (true) {
            test_count++;
            if (test_count % 45000 == 0) {
                total_minutes += 5;
                cerr << "Done " << test_count << '\n';
                cerr << total_minutes << " minutes passed\n";
            }
            clock_t CURRENT_TIME = clock();
            if (CURRENT_TIME - START_TIME > Time * CLOCKS_PER_SEC)
                break;
            const int UPPER = 25;
            auto gen_pair = [&]() -> pair<int, int> {
                //return {rnd() % UPPER - UPPER / 2, rnd() % UPPER - UPPER / 2};
                return {rnd() % UPPER, rnd() % UPPER};
            };
            // adjacency matrix
            // todo remove weights
            vector<vector<pair<int, int>>> weights(n, vector<pair<int, int>>(n));

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    weights[i][j] = gen_pair();
                }
            vector<int> best[2];
            vector<vector<pair<int, int>>> mat;
            int N = strategy[0].size(), M = strategy[1].size();

            mat.resize(N, vector<pair<int, int>>(M));
            for (int t = 0; t < 2; t++)
                best[t].resize(strategy[t].size(), numeric_limits<int>::max());
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++) {
                    auto x = get_score(weights, i, j);
                    mat[i][j] = x;
                    if (best[0][i] > x.second)
                        best[0][i] = x.second;
                    if (best[1][j] > x.first)
                        best[1][j] = x.first;
                }
            pair<int, int> nash = {-1, -1};
            for (int i = 0; i < N && nash.first == -1; i++)
                for (int j = 0; j < M; j++) {
                    auto x = mat[i][j];
                    if (x.second == best[0][i] && x.first == best[1][j]) {
                        nash = {i, j};
                        break;
                    }
                }
            if (nash.first == -1) {
                cerr << "FOUND " << test_count << "\n";
                for (int i = 0; i < n; i++)
                    for (int j : edges[i])
                        cout << i + 1 << ' ' << j + 1 << ' ' << weights[i][j].first << ' ' << weights[i][j].second << '\n';
                return;
            }
        }
        cerr << "Ran " << test_count << " tests\n";
    }
};
