#pragma once

#include "Solver.hpp"
#include <bits/stdc++.h>
using namespace std;

// Tries assigning random integers and checks for Nash Equilibrium
struct LinearSolver : Solver {

    LinearSolver(atomic<int> &X_, atomic_bool &require_exit, int thread_id) : Solver(X_, require_exit, thread_id) {}

    int get_score(vector<int> const &str, int st = -1) {

        vector<int> out(n, -1);
        for (int t = 0; t < players; t++) {
            int id = str[t];
            for (int i = 0; i < (int)assigned[t].size(); i++)
                out[assigned[t][i]] = strategy[t][id][i];
        }
//        cerr << "Current strategy:\n";
//        for (int i = 0; i < n; i++)
//            cerr << out[i] << ' ';
//        cerr << '\n';
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
            return terminal;
        }
        assert(term_id[v] != -1);
        return term_id[v];
    }

    int test_count = 0;

    bool found_nash;

    bool found_non_nash;

    bool preprocess(bool print) {
        found_nash = false;
        found_non_nash = false;

        int total = total_encode(-1);

        int left = 0, right = n;

        bool uniform_nash = false;
        assert(!uniform_nash);

        if (!uniform_nash) {
            left = start;
            right = start + 1;
        }

        constraints.clear();

        for (int st = left; st < right; st++)
            if (!edges[st].empty()) {
                /*vector<vector<int>> best(players);
                for (int i = 0; i < players; i++) {
                    best[i].resize(total_encode(i), numeric_limits<int>::min());
                }*/

                vector<int> mat(total);

                for (int i = 0; i < total; i++) {
                    auto str = decode(i, -1);
                    auto x = get_score(str, st);
                    mat[i] = x;
                    //cerr << "strat " << i << ": " << x << '\n';
                    /*for (int j = 0; j < players; j++) {
                        int id = encode(str, j);
                        if (best[j][id] < x[j]) {
                            best[j][id] = x[j];
                        }
                    }*/
                }
                for (int i = 0; i < total; i++) {
                    auto str = decode(i, -1);
                    constraints.push_back({mat[i], {}});
                    vector<pair<int, int>> who;
                    for (int j = 0; j < players; j++) {
                        for (int k = 0; k < strategy[j].size(); k++) {
                            int old_value = str[j];
                            str[j] = k;
                            who.push_back({j, get_score(str, st)});
                            str[j] = old_value;
                        }
                    }
                    //cerr << "strategy " << i << ": " << mat[i] << ":\n";
                    //for (auto [pl, j] : who)
                    //    cerr << "(" << pl << ", " << j << ") ";
                    //cerr << '\n';
                    sort(who.begin(), who.end());
                    who.resize(unique(who.begin(), who.end()) - who.begin());
                    if (bad_cycles) {
                        vector<pair<int, int>> who_copy;
                        for (auto j : who)
                            if (j.second != terminal)
                                who_copy.push_back(j);
                        swap(who, who_copy);
                    }
                    {
                        vector<pair<int, int>> who_copy;
                        for (auto j : who)
                            if (j.second != mat[i])
                                who_copy.push_back(j);
                        swap(who, who_copy);
                    }
                    if (who.empty()) {
                        found_nash = true;
                    }
                    constraints.back().second = who;
                    if (bad_cycles && constraints.back().first == terminal)
                        constraints.pop_back();
                }
            }
        bool require_bad_cycle = false;
        bad_mask = 0;

        PLAYERS = 0;
        if (require_bad_cycle) {
            for (int i = term_cyc_seg.first; i < term_cyc_seg.second; i++) {
                assert(term_id[i] != -1);
                bad_mask ^= (1 << term_id[i]);
                for (int player = 0; player < PLAYERS; player++)
                    constraints.emplace_back(terminal, vector{make_pair(player, term_id[i])});
            }
        }
        if (found_nash) {
            if (print)
                cerr << "This example has nash equilibrium because some strategy can't be exited\n";
            return false;
        }
        sort(constraints.begin(), constraints.end(), [&](auto x, auto y) {
            return make_pair(x.first, x.second.size()) < make_pair(y.first, y.second.size());
        });

        {
            vector<pair<int, vector<pair<int, int>>>> c_copy;
            for (auto[node, con] : constraints) {
                bool good = 1;
                for (auto[_, c1] : c_copy)
                    if (_ == node) {
                        set<pair<int, int>> g1{c1.begin(), c1.end()};
                        bool inside = 1;
                        for (auto i : con)
                            if (!g1.count(i)) {
                                inside = 0;
                                break;
                            }
                        if (inside) {
                            good = 0;
                            break;
                        }
                    }
                if (good) {
                    c_copy.emplace_back(node, con);
                }
            }
            swap(constraints, c_copy);
        }

        if (require_bad_cycle) {
            vector<pair<int, vector<pair<int, int>>>> c_copy;
            for (auto [node, con] : constraints)
                if (node == terminal && con.size() == 1) {
                    int v = ob_term[con[0].second];
                    if (con[0].first < PLAYERS && v >= term_cyc_seg.first && v < term_cyc_seg.second);
                    else c_copy.emplace_back(node, con);
                } else {
                    c_copy.emplace_back(node, con);
                }
            swap(constraints, c_copy);
        }

        {
            edge_to_constraint.resize(players);
            for (auto &i : edge_to_constraint) {
                i.resize(terminal + 1);
                for (auto &j : i) {
                    j.clear();
                    j.resize(terminal + 1);
                    for (auto z : j)
                        assert(z.empty());
                }
            }

            back_edge_to_constraint.resize(players);
            for (auto &i : back_edge_to_constraint) {
                i.resize(terminal + 1);
                for (auto &j : i) {
                    j.clear();
                    j.resize(terminal + 1);
                    for (auto z : j)
                        assert(z.empty());
                }
            }

            int index = 0;
            for (auto [t_0, pairs] : constraints) {
                for (auto [p_i, t_i] : pairs) {
                    edge_to_constraint[p_i][t_0][t_i].push_back(index);
                    back_edge_to_constraint[p_i][t_i][t_0].push_back(index);
                }
                index++;
            }
        }

        weights_order.clear();
        weights_order.resize(players, vector<int>(terminal + 1));
        unsatisfied_constraints.resize(constraints.size());
        for (int i = 0; i < constraints.size(); i++) {
            unsatisfied_constraints[i] = constraints[i].second.size();
        }
        satisfied_constraints.clear();
        satisfied_constraints.resize(constraints.size(), 0);

        present_in_mask.resize((1 << (terminal + 1)));
        for (int mask = 0; mask < present_in_mask.size(); mask++) {
            present_in_mask[mask].clear();
            for (int i = 0; i <= terminal; i++)
                if ((1 << i) & mask)
                    present_in_mask[mask].push_back(i);
        }

        players_order.resize(players);
        iota(players_order.begin(), players_order.end(), 0);
        {
            vector<int> counter(players);
            for (int i = 0; i < players; i++)
                for (int x = 0; x <= terminal; x++)
                    for (int y = 0; y <= terminal; y++)
                        counter[i] += edge_to_constraint[i][x][y].size();
            sort(players_order.begin(), players_order.end(), [&](int i, int j) {
                return counter[i] > counter[j];
            });
        }

        already_satisfied.resize(players);
        for (auto &i : already_satisfied) {
            i.resize(terminal + 1);
            for (auto &x : i)
                x.clear();
        }

        count_bad_cycles = 0;

        return true;
    }

    int PLAYERS;

    bool bad_cycles = false;
    const int max_bad_cycle = 2, max_bad_cycle_id = 3;

    //<t_0, vector<p_i, t_j>> - t_0 is less than at least one of <p_i, t_j>
    vector<pair<int, vector<pair<int, int>>>> constraints;

    // array<array<array<vector<int>, terminal+1>, terminal+1>, players>
    // <p, v, u> = z means constraint z has v < u for player p
    vector<vector<vector<vector<int>>>> edge_to_constraint, back_edge_to_constraint;

    // (players, terminal + 1)
    vector<vector<int>> weights_order;

    vector<int> satisfied_constraints;

    vector<int> unsatisfied_constraints;

    vector<vector<int>> present_in_mask;

    vector<int> players_order;

    vector<vector<vector<int>>> already_satisfied;

    int count_bad_cycles;

    int bad_mask;

    void rec_brute(int player_id, int id, int left_mask) {
        if (player_id == players) {

            vector<vector<int>> weights(terminal + 1, vector<int>(players));
            for (int i = 0; i < players; i++)
                for (int j = 0; j <= terminal; j++)
                    weights[weights_order[i][j]][i] = terminal - j;
            bool good_thing = 1;
            if (false) {
                auto q = weights.back();
                sort(q.begin(), q.end());
                vector<int> who = {2, 4};
                for (int i = 0; i < who.size(); i++)
                    if (q[q.size() - 1 - i] > who[who.size() - 1 - i]) {
                        good_thing = 0;
                        break;
                    }
            }
            if (good_thing) {
                found_non_nash = true;
                cerr << "FOUND\n";
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
            }
            return;
        }

        // todo change player_id to player?
        int player = players_order[player_id];
        if (present_in_mask[left_mask].size() == 1) {
            weights_order[player][id] = present_in_mask[left_mask][0];
            rec_brute(player_id + 1, 0, (1 << (terminal + 1)) - 1);
            return;
        }
        auto try_value = [&](int x, bool is_great_value) {
            bool bad = 0;
            for (int y : present_in_mask[left_mask]) {
                if (!is_great_value) {
                    for (int i : edge_to_constraint[player][x][y])
                        if (!satisfied_constraints[i]) {
                            unsatisfied_constraints[i]--;
                            if (unsatisfied_constraints[i] == 0)
                                bad = 1;
                        }
                }
                for (int i : back_edge_to_constraint[player][x][y])
                    if (!satisfied_constraints[i]) {
                        satisfied_constraints[i] = 1;
                        already_satisfied[player_id][id].push_back(i);
                    }
            }

            if (!bad) {
                weights_order[player][id] = x;
                rec_brute(player_id, id + 1, left_mask ^ (1 << x));
            }

            for (int i : already_satisfied[player_id][id])
                satisfied_constraints[i] = 0;
            already_satisfied[player_id][id].clear();

            if (!is_great_value) {
                for (int y : present_in_mask[left_mask])
                    for (int i : edge_to_constraint[player][x][y])
                        if (!satisfied_constraints[i]) {
                            unsatisfied_constraints[i]++;
                        }
            }
        };

        int great_value = -1;

        for (int x : present_in_mask[left_mask])
            //if (x != terminal || (player >= PLAYERS && !bad_cycles && count_bad_cycles > players)) {
            if (x != terminal) {
                bool bad = 0;
                for (int y : present_in_mask[left_mask]) {
                    for (int i : edge_to_constraint[player][x][y])
                        if (!satisfied_constraints[i]) {
                            bad = 1;
                            break;
                        }
                    if (bad)
                        break;
                }
                if (!bad) {
                    great_value = x;
                    break;
                }
            }
        if (great_value != -1) {
            try_value(great_value, true);
        } else {
            for (int x : present_in_mask[left_mask])
                if (x != terminal ||
                //(player >= PLAYERS && !bad_cycles && count_bad_cycles < max_bad_cycle && present_in_mask[left_mask].size() - 1 <= max_bad_cycle_id)) {
                (player >= PLAYERS && !bad_cycles && count_bad_cycles < max_bad_cycle &&
                        present_in_mask[left_mask].size() - 1 <= max_bad_cycle_id)) {
                //|| (player == 0 && count_bad_cycles < max_bad_cycle && present_in_mask[left_mask].size() - 1 <= 2)) {
                //(!bad_cycles && (player >= PLAYERS || (left_mask & bad_mask) == 0))) {
                    if (x == terminal)
                        count_bad_cycles++;
                    try_value(x, false);
                    if (found_non_nash)
                        return;
                    if (x == terminal)
                        count_bad_cycles--;
                }
        }
        //todo dont allocate memory
    }

    // Time - number of seconds to run
    bool solve(bool print = true, int Time = 600000, int random_seed = -1) {
        if (random_seed == -1)
            random_seed = chrono::high_resolution_clock::now().time_since_epoch().count();
        mt19937 rnd(random_seed);

        clock_t START_TIME = clock();
        {
            long long fact = 1;
            for (int i = 2; i <= terminal + !bad_cycles; i++)
                fact *= i;
            long long res = 1;
            for (int i = 0; i < players; i++)
                res *= fact;
            if (print)
                cerr << "total of " << res << " orders\n";
        }

        if (false) {
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
        if (false) {
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

        if (!preprocess(print))
            return false;

//        for (auto [x, y] : constraints) {
//            cerr << x << ' ';
//            for (auto [z, o] : y)
//                cerr << "(" << z << ", " << o << ") ";
//            cerr << '\n';
//        }

        rec_brute(0, 0, (1 << (terminal + 1)) - 1);

        return found_non_nash;
    }
};
