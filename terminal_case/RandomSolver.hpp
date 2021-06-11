#pragma once

#include "Solver.hpp"
#include <bits/stdc++.h>
using namespace std;

// Tries assigning random integers and checks for Nash Equilibrium
struct RandomSolver : Solver {

    RandomSolver(atomic<int> &X_, atomic_bool &require_exit, int thread_id) : Solver(X_, require_exit, thread_id) {}

    vector<int> get_score(vector<vector<int>> const &weights, vector<int> const &str, int st = -1) {

        vector<int> out(n, -1);
        for (int t = 0; t < players; t++) {
            int id = str[t];
            for (int i = 0; i < (int)assigned[t].size(); i++)
                out[assigned[t][i]] = strategy[t][id][i];
        }
        vector<bool> visited(n);
        int v = st;
        if (v == -1)
            v = start;
        while (!visited[v] && !edges[v].empty()) {
            visited[v] = true;
            v = out[v];
            assert(v >= 0);
        }
        if (visited[v]) {
            return weights[terminal];
        }
        assert(term_id[v] != -1);
        return weights[term_id[v]];
    }

    int test_count = 0;

    bool found_nash;

    int total_good;

    void solve(vector<vector<int>> weights) {
        test_count++;
        if (test_count % 50000 == 0) {
            cerr << "Done " << test_count << '\n';
            //cerr << total_minutes << " minutes passed\n";
        }
        int total = total_encode(-1);

        vector<int> nash_eq(total);

        int total_times = 0;

        found_nash = false;

        int left = 0, right = n;

        bool uniform_nash = true;

        if (!uniform_nash) {
            left = start;
            right = start + 1;
        }

        for (int st = left; st < right; st++)
            if (!edges[st].empty()) {
                vector<vector<int>> best(players);
                for (int i = 0; i < players; i++) {
                    best[i].resize(total_encode(i), numeric_limits<int>::min());
                }

                vector<vector<int>> mat(total);

                for (int i = 0; i < total; i++) {
                    auto str = decode(i, -1);
                    auto x = get_score(weights, str, st);
                    mat[i] = x;
                    for (int j = 0; j < players; j++) {
                        int id = encode(str, j);
                        if (best[j][id] < x[j]) {
                            best[j][id] = x[j];
                        }
                    }
                }
                found_nash = false;
                for (int i = 0; i < total; i++)
                    if (nash_eq[i] == total_times) {
                        auto str = decode(i, -1);
                        auto &x = mat[i];
                        bool good = true;
                        for (int j = 0; j < players; j++) {
                            int id = encode(str, j);
                            if (best[j][id] != x[j]) {
                                assert(best[j][id] > x[j]);
                                good = false;
                                break;
                            }
                        }
                        if (good) {
                            found_nash = true;
                            nash_eq[i]++;
                            //break;
                        }
                    }

                total_times++;
                if (!found_nash) {
                    break;
                }
            }

        if (!found_nash) {
            total_good++;
/*#define cerr cout
            cerr << "FOUND " << test_count << "\n";
            cerr << "graph:\n";
            cerr << n << ' ' << players << ' ' << start + 1 << '\n';
            for (int i = 0; i < n; i++) {
                cerr << assignment[i] + 1 << ' ';
            }
            cerr << '\n';
            for (int i = 0; i < n; i++) {
                for (int j : edges[i]) {
                    cerr << i + 1 << ' ' << j + 1 << '\n';
                }
            }
            cerr << "\nexample:\n";
            for (int i = 0; i <= terminal; i++) {
                if (i < terminal) {
                    cerr << "node " << ob_term[i] + 1 << '\n';
                } else {
                    cerr << "cycle weights\n";
                }
                for (int j : weights[i]) {
                    cerr << j << ' ';
                }
                cerr << '\n';
            }
            cerr << '\n';
#undef cerr*/
            cerr << "FOUND " << test_count << "\n";
            cerr << "graph:\n";
            cerr << n << ' ' << players << ' ' << start + 1 << '\n';
            for (int i = 0; i < n; i++) {
                cerr << assignment[i] + 1 << ' ';
            }
            cerr << '\n';
            for (int i = 0; i < n; i++) {
                for (int j : edges[i]) {
                    cerr << i + 1 << ' ' << j + 1 << '\n';
                }
            }
            cerr << "\nexample:\n";
            for (int i = 0; i <= terminal; i++) {
                if (i < terminal) {
                    cerr << "node " << ob_term[i] + 1 << '\n';
                } else {
                    cerr << "cycle weights\n";
                }
                for (int j : weights[i]) {
                    cerr << j << ' ';
                }
                cerr << '\n';
            }
            cerr << '\n';
            //exit(0);
        }
    }

    const vector<int> cycle_ids = {0};

    bool bad_cycles = false;

    vector<int> cur_cycles;

    void rec_brute(vector<vector<int>> &weights, int id) {
        if (id == players) {
            solve(weights);
            return;
        }
        vector<int> p(terminal + 1);
        iota(p.begin(), p.end(), 0);
        reverse(p.begin(), p.end());
        do {
            bool allowed = 1;
            int to_rem = -1;
            if (p.back() != 0) {
                if (cur_cycles.size() == cycle_ids.size()) {
                    allowed = 0;
                } else {
                    to_rem = lower_bound(cur_cycles.begin(), cur_cycles.end(), p.back()) - cur_cycles.begin();
                    cur_cycles.insert(cur_cycles.begin() + to_rem, p.back());
                    for (int i = 0; i + 1 < cur_cycles.size(); i++)
                        assert(cur_cycles[i] <= cur_cycles[i + 1]);
                    for (int j = 0; j < cur_cycles.size(); j++)
                        if (cur_cycles[cur_cycles.size() - 1 - j] > cycle_ids[cycle_ids.size() - 1 - j]) {
                            allowed = 0;
                            break;
                        }
                }
            }
            if (allowed) {
                for (int i = 0; i <= terminal; i++)
                    weights[i][id] = p[i];
                rec_brute(weights, id + 1);
                //if(!found_nash)
                //    return;
            }
            if (to_rem != -1)
                cur_cycles.erase(cur_cycles.begin() + to_rem);
        } while (prev_permutation(p.begin(), p.end()));
    }

    // Time - number of seconds to run
    bool solve(bool print = true, int Time = 600000, int random_seed = -1) {
        if (random_seed == -1)
            random_seed = chrono::high_resolution_clock::now().time_since_epoch().count();
        mt19937 rnd(random_seed);

        clock_t START_TIME = clock();
        test_count = 0;
        {
            long long fact = 1;
            for (int i = 2; i <= terminal + !bad_cycles; i++)
                fact *= i;
            long long res = 1;
            for (int i = 0; i < players; i++)
                res *= fact;
            if (true or print)
                cerr << "total of " << res << " orders\n";
        }

        {
            vector<bool> ex(players);
            for (int i : assignment)
                if (i != -1)
                    ex[i] = 1;
            for (bool i : ex)
                if (!i) {
                    if (print) {
                        cerr << "Some player is not present\n";
                    }
                    return false;
                }
        }
        {
            for (int i = 0; i < n; i++)
                if (!edges[i].empty() && edges[i].size() == 1) {
                    if (print) {
                        cerr << "Node " << to_string(i) << " has only 1 edge\n";
                    }
                    return false;
                }
        }
        if (precheck)
            return true;

        vector<vector<int>> weights(terminal + 1, vector<int>(players));
        bool random = true;
        if (random) {
            total_good = 0;
            for (int i = 0; i + 1 < cycle_ids.size(); i++)
                assert(cycle_ids[i] <= cycle_ids[i + 1]);
            cur_cycles.clear();
            if (bad_cycles)
                assert(cycle_ids.empty());
            rec_brute(weights, 0);
            if (print) {
                cerr << "deterministic not found\n";
            }
            cerr << "Total good is " << total_good << '\n';
        } else {
            while (true) {
                clock_t CURRENT_TIME = clock();
                if (CURRENT_TIME - START_TIME > Time * CLOCKS_PER_SEC)
                    break;
                const int UPPER = 1e9;
                auto gen_weights = [&]() -> vector<int> {
                    vector<int> res(players);
                    for (int &i : res)
                        i = rnd() % UPPER;

                    return res;
                };


                // terminal - weights of a cycle
                for (int i = 0; i <= terminal; i++) {
                    weights[i] = gen_weights();
                }
                //
                if (bad_cycles) {
                    for (int i = 0; i < players; i++)
                        weights[terminal][i] = -1;
                }
                solve(weights);
            }
            cerr << "random not found\n";
        }
        if (print)
            cerr << "Ran " << test_count << " tests\n";
        return false;
    }
};
