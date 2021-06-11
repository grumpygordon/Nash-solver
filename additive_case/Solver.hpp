#pragma once

#include <vector>
#include <iostream>
using namespace std;

// call solve() of any Solver successor
// Currently either RandomSolver or SlowSolver
struct Solver {
    // n - number of nodes
    // edges - adjacency list
    // start - source
    // finish - sink

    Solver() = default;

    virtual void solve() = 0;

    int m;
    vector<vector<int>> edges;
    vector<vector<int>> edges_id;

    vector<int> assignment;

    int start, finish, n;

    // reads input from stdin
    // for exact format call help()

    // list of assigned nodes to a player
    vector<int> assigned[2];

    // list of assigned edges ends[u from (v->u) pair] for each node which is assigned to a player
    // for each strategy
    // same order as assigned
    vector<vector<int>> strategy[2];
    vector<vector<int>> strategy_id[2];

private:
    void brute_strats(vector<int> &cur, vector<int> &cur_id, int t, int id) {
        if (id == (int)assigned[t].size()) {
            strategy[t].push_back(cur);
            strategy_id[t].push_back(cur_id);
            return;
        }
        cur.push_back(-1);
        cur_id.push_back(-1);
        bool tried = false;
        int v = assigned[t][id];
        for (int i = 0; i < (int)edges[v].size(); i++) {
            tried = true;
            cur.back() = edges[v][i];
            cur_id.back() = i;
            brute_strats(cur, cur_id, t, id + 1);
        }
        if (!tried)
            brute_strats(cur, cur_id, t, id + 1);
        cur.pop_back();
        cur_id.pop_back();
    }

    bool make_checks() {
        assert(n != 0);
        {
            // checking that number of strategies for one player is lower than BOUND = 10^3
            int N = 1, M = 1;
            const int BOUND = 1e3;

            for (int i = 0; i < n; i++)
                if (assignment[i] == 0) {
                    N *= edges[i].size();
                    if (N > BOUND)
                        return false;
                } else {
                    assert(assignment[i] == 1);
                    M *= edges[i].size();
                    if (M > BOUND)
                        return false;
                }
            return true;
            return N * M <= 150;
        }
    }

    vector<bool> visited;

    void dfs(int v) {
        visited[v] = 1;
        for (int i : edges[v])
            if (!visited[i])
                dfs(i);
    }

    bool precalc() {
        visited.assign(n, 0);
        dfs(start);
        if (!visited[finish])
            return 0;
        // nodes assigned for first and second players respectively
        for (int t = 0; t < 2; t++) {
            assigned[t].clear();
            strategy[t].clear();
            strategy_id[t].clear();
        }

        {
            for (int i = 0; i < n; i++) {
                int t = assignment[i];
                assigned[t].push_back(i);
            }
        }
        edges_id.assign(edges.size(), {});
        m = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < (int)edges[i].size(); j++) {
                edges_id[i].push_back(m);
                m++;
            }
        }
        for (int t = 0; t < 2; t++) {
            vector<int> cur;
            vector<int> cur_id;
            cur.reserve(assigned[t].size());
            cur_id.reserve(assigned[t].size());
            brute_strats(cur, cur_id, t, 0);
        }
        cerr << "We have " << strategy[0].size() << ' ' << strategy[1].size() << " strategies\n";
        cerr << "Total of " << strategy[0].size() * (long long)strategy[1].size() << " strategies\n";
        return strategy[0].size() * strategy[1].size() <= 200;
    }

public:

    void gen_tree() {
        ifstream cin("../a.in");
        int prefix, cycle_len;
        cin >> prefix >> cycle_len;
        //if (thread_id == 0)
            cerr << "n c " << prefix << ' ' << cycle_len << '\n';
        vector<int> cyc_ass(cycle_len);
        for (int &i : cyc_ass) {
            cin >> i;
            i--;
        }
        vector<pair<int, int>> cyc_edges;
        int cycle_nodes = cycle_len;
        {
            int sz;
            cin >> sz;
            cyc_edges.resize(sz);
            for (auto &[v, u] : cyc_edges) {
                cin >> v >> u;
                v--;
                u--;
            }
        }
        n = prefix + cycle_nodes;
        start = 0;
        finish = n - 1;
        assignment.resize(n, -1);
        int IT = 0;
        //auto perm_n = generate_permutations(n);
        //auto perm_c = generate_permutations(cycle_len);
        edges.resize(n);
        mt19937 rnd(2283221337 * 1 + int('p' + 'l' + 'z' + 'h' + 'e' + 'l' + 'p'));
        while (true) {
            edges.assign(n, {});
            vector<int> parent(prefix);
            for (int i = 1; i < prefix; i++) {
                parent[i] = rnd() % i;
                edges[parent[i]].push_back(i);
            }
            for (int i = 0; i < cycle_len; i++) {
                assignment[prefix + i] = cyc_ass[i];
                //assignment[prefix + i] = rnd() % (players - 1);
            }
            for (auto [v, u] : cyc_edges)
                edges[prefix + v].push_back(prefix + u);
            int skoka = 0;
            for (int i = 0; i < prefix; i++) {
                for (int j = i + 1; j < prefix; j++)
                    if (parent[j] != i) {
                        if (rnd() % (prefix + 2) == 0)
                            edges[i].push_back(j);
                    }
                for (int j = 0; j < cycle_len - 1; j++) {
                    if (rnd() % ((prefix + 1)) == 0) {
                        edges[i].push_back(prefix + j);
                        skoka++;
                    }
                }
            }
            int additional = 0;
            for (int i = 0; i < prefix; i++) {
                if (edges[i].empty()) {
                    edges[i].push_back(i + 1 + rnd() % (prefix + cycle_len - 1 - i));
                }
                if (additional < 2 && rnd() % max(2, prefix + 1) == 0) {
                    edges[i].push_back(finish);
                }
                assignment[i] = rnd() % 2;
            }
            if (!make_checks())
                continue;
            if (!precalc())
                continue;
            solve();
            if (IT % 50 == 0)
                cerr << "Done " + to_string(IT) + "\n";
            IT++;
            // never happens
            if (IT == -1)
                break;
        }
        //cerr << "There is " << good_count << " non nash equilibriums in total\n";
        //cerr << "ratio is " << good_count / (double)(min(IT, r_start) - max(0, l_start)) << '\n';
    }


    void gen_tree_cycle() {
        ifstream cin("../a.in");
        int prefix, cycle_len;
        cin >> prefix >> cycle_len;
            cerr << "n p s " << prefix << ' ' << cycle_len << '\n';
        start = 0;
        finish = prefix + cycle_len;
        int IT = 0;
        //auto perm_n = generate_permutations(n);
        //auto perm_c = generate_permutations(cycle_len);
        mt19937 rnd(12540 * 1 + 223);
        while (true) {
            n = prefix + cycle_len + 1;
            assignment.assign(n, 0);
            edges.assign(n, {});
            vector<int> parent(prefix);
            for (int i = 1; i < prefix; i++) {
                parent[i] = rnd() % i;
                edges[parent[i]].push_back(i);
            }
            for (int i = 0; i < cycle_len; i++) {
                {
                    int o = rnd() % (cycle_len - 1);
                    if (o < i)
                        edges[prefix + i].push_back(prefix + o);
                    else
                        edges[prefix + i].push_back(prefix + o + 1);
                }
                for (int j = 0; j < cycle_len; j++) {
                    // todo change that
                    if (i != j && j != edges[prefix + i][0] && rnd() % (cycle_len) == 0) {
                        edges[prefix + i].push_back(prefix + j);
                    }
                }
                assignment[prefix + i] = rnd() % 2;
                if (rnd() % ((cycle_len + 1)) == 0)
                    edges[prefix + i].push_back(prefix + cycle_len);
            }
            for (int i = 0; i < prefix; i++) {
                for (int j = i + 1; j < prefix; j++)
                    if (parent[j] != i) {
                        if (rnd() % (prefix + 2) == 0)
                            edges[i].push_back(j);
                    }
                for (int j = 0; j < cycle_len; j++) {
                    if (rnd() % (cycle_len + 1) == 0)
                        edges[i].push_back(prefix + j);
                }
            }
            int additional = 0;
            for (int i = prefix - 1; i >= 0; i--) {
                if (edges[i].empty()) {
                    edges[i].push_back(i + 1 + rnd() % (prefix + cycle_len - 1 - i));
                }
                if (i > 0 && additional < 2 && rnd() % max(1, prefix + 1) == 0) {
                    edges[i].push_back(n - 1);
                }
                assignment[i] = rnd() % 2;
            }
            make_checks();
            if (!precalc())
                continue;
            solve();
            if (IT % 50 == 0)
                cerr << "Done " + to_string(IT)+ "\n";
            IT++;
            if (IT == -1)
                break;
        }
        //cerr << "There is " << good_count << " non nash equilibriums in total\n";
        //cerr << "ratio is " << good_count / (double)(min(IT, r_start) - max(0, l_start)) << '\n';
    }

    void read() {
        cin >> n;
        cin >> start >> finish;
        start--;
        finish--;
        edges.resize(n);
        assignment.resize(n);
        for (int i = 0; i < n; i++) {
            cin >> assignment[i];
            assert(assignment[i] >= 1 && assignment[i] <= 2);
            assignment[i]--;
        }
        fill(edges.begin(), edges.end(), vector<int>());
        int type;
        cin >> type;
        if (type == 1) {
            int edges_number;
            cin >> edges_number;
            for (int i = 0; i < edges_number; i++) {
                int v, u;
                cin >> v >> u;
                assert(min(v, u) >= 1 && max(v, u) <= n);
                edges[v - 1].push_back(u - 1);
            }
        } else {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    int w;
                    cin >> w;
                    if (w == 1)
                        edges[i].push_back(j);
                }
        }
        make_checks();
        precalc();
    }

    static void help() {
        cerr << "# Enter following:\n";
        cerr << "n\n";
        cerr << "s t\n";
        cerr << "# each a_i in {1, 2}\n";
        cerr << "a_1 a_2 ... a_n\n";
        cerr << "# Now enter 'type': how you would like to enter edges\n";
        cerr << "# if type = 1, then adjacency list in the form:\n";
        cerr << "m\n";
        cerr << "v_1 u_1\n";
        cerr << "v_2 u_2\n";
        cerr << "...\n";
        cerr << "v_m u_m\n";
        cerr << "# if type = 2, then adjacency matrix in the form,\n";
        cerr << "# where 1 corresponds to existing edge:\n";
        cerr << "e_1_1 e_1_2 ... e_1_n\n";
        cerr << "e_2_1 e_2_2 ... e_2_n\n";
        cerr << "...   ...   ...   ...\n";
        cerr << "e_n_1 e_n_2 ... e_n_n\n";
        cerr << "\nExample:\n";
        cerr << "4\n"
                "1 4\n"
                "1 2 2 1\n"
                "\n"
                "1\n"
                "4\n"
                "1 2\n"
                "1 3\n"
                "2 4\n"
                "3 4\n";
        cerr << "\nSame example:\n";
        cerr << "4\n"
                "1 4\n"
                "1 2 2 1\n"
                "\n"
                "2\n"
                "0 1 1 0\n"
                "0 0 0 1\n"
                "0 0 0 1\n"
                "0 0 0 0\n";

    }
};
