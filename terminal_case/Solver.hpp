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
    // players - number of players

    virtual bool solve(bool print = true, int = 0, int = 0) = 0;

    int m;
    // (m)
    vector<vector<int>> edges;

    // id of player, whom node is assigned to
    // -1 in case of terminal node
    // (n)
    vector<int> assignment;

    // term_id[v] = x means node v is terminal and has internal index x
    // -1 for non terminal node
    // (n)
    vector<int> term_id;

    // index of node by index of term_id
    vector<int> ob_term;

    // number of terminal nodes
    int terminal;

    int start, n;
    int players;

    // reads input from stdin
    // for exact format call help()

    // list of assigned nodes to a player
    // (players)
    vector<vector<int>> assigned;

    // list of assigned edges ends[u from (v->u) pair] for each node which is assigned to a player
    // for each strategy
    // same order as assigned
    // (players)
    vector<vector<vector<int>>> strategy;

private:

    int make_checks() {
        assert(n != 0);
        {
            //checking that terminal node is unassigned
            for (int i = 0; i < n; i++) {
                bool a = (assignment[i] == -1);
                bool b = (edges[i].empty());
                if (a != b) {
                    cerr << "node " << i + 1 << " is bad\n";
                }
                assert(a == b);
            }
        }
        bool failed = 0;
        int TOT = 1;
        {
            // checking that number of strategies for one player is lower than BOUND = 10^3
            vector<int> str(players, 1);
            const int BOUND = 1e6;

            for (int i = 0; i < n; i++) {
                int id = assignment[i];
                if (id != -1) {
                    TOT *= edges[i].size();
                    //str[id] *= edges[i].size();
                    if (TOT > BOUND)
                        failed = 1;
                    //assert(TOT <= BOUND);
                }
            }
        }
        if (failed)
            return TOT;
        return 0;
    }

    void brute_strats(vector<int> &cur, int t, int id) {
        if (id == (int)assigned[t].size()) {
            strategy[t].push_back(cur);
            return;
        }
        cur.push_back(-1);
        bool tried = false;
        int v = assigned[t][id];
        for (int u : edges[v]) {
            tried = true;
            cur.back() = u;
            brute_strats(cur, t, id + 1);
        }
        assert(tried);
        cur.pop_back();
    }

    void precalc(bool print = true) {
        term_id.resize(n);
        strategy.clear();
        assigned.clear();
        strategy.resize(players);
        assigned.resize(players);
        {
            terminal = 0;
            ob_term.clear();
            for (int i = 0; i < n; i++) {
                int t = assignment[i];
                if (t == -1) {
                    term_id[i] = terminal;
                    ob_term.push_back(i);
                    terminal++;
                } else {
                    term_id[i] = -1;
                    assigned[t].push_back(i);
                }
            }
        }
        for (int t = 0; t < players; t++) {
            vector<int> cur;
            cur.reserve(assigned[t].size());
            brute_strats(cur, t, 0);
            assert(!strategy[t].empty());
        }
        if (print) {
            long long total = 1;
            cerr << "We have";
            for (auto &i : strategy) {
                cerr << ' ' << i.size();
                total *= i.size();
            }
            cerr << " strategies\n";
            cerr << "Total of " << total << " strategies\n";
            assert(total <= int(1e6));
        }
    }

public:

    // assign - 0 if we should decide, otherwise id of player
    // e1 - must have edges
    // e2 - edge space to brute

    int counter = 0;

    bool precheck = false;

    atomic<int> &X;

    atomic_bool &require_exit;

    int thread_id;

    Solver() = default;

    Solver(atomic<int> &T, atomic_bool &_require_exit, int y) : X(T), require_exit(_require_exit), thread_id(y) {}

    vector<vector<vector<int>>> generate_permutations(int n) {
        vector<vector<vector<int>>> res(n + 1);
        for (int mask = 0; mask < (1 << n); mask++) {
            vector<int> o;
            for (int i = 0; i < n; i++)
                if ((1 << i) & mask)
                    o.push_back(i);
            res[o.size()].push_back(o);
        }
        return res;
    }

    void gen_tree() {
        ifstream cin("../a.in");
        int prefix, cycle_len;
        cin >> prefix >> cycle_len >> players;
        if (thread_id == 0)
            cerr << "n p s " << prefix << ' ' << cycle_len << ' ' << players << '\n';
        vector<int> cyc_ass(cycle_len);
        for (int &i : cyc_ass) {
            cin >> i;
            i--;
        }
        vector<pair<int, int>> cyc_edges;
        int cycle_nodes = 0;
        {
            int sz;
            cin >> sz;
            cyc_edges.resize(sz);
            for (auto &[v, u] : cyc_edges) {
                cin >> v >> u;
                v--;
                u--;
            }
            //for (int i = 0; i < cycle_len; i++)
            //    cyc_edges.push_back({i, (i + 1) % cycle_len});
            vector<int> different;
            for (auto [v, u] : cyc_edges) {
                different.push_back(v);
                different.push_back(u);
            }
            sort(different.begin(), different.end());
            different.resize(unique(different.begin(), different.end()) - different.begin());
            map<int, int> diff;
            for (int i = 0; i < different.size(); i++)
                diff[different[i]] = i;
            cycle_nodes = different.size();
            for (auto &[v, u] : cyc_edges) {
                v = diff[v];
                u = diff[u];
            }
        }
        start = 0;
        n = prefix + cycle_nodes;
        term_cyc_seg = {prefix + cycle_len, n};
        //term_cyc_seg.first += 100;
        assignment.resize(n, -1);
        int IT = 0;
        //auto perm_n = generate_permutations(n);
        //auto perm_c = generate_permutations(cycle_len);
        edges.resize(n);
        mt19937 rnd(2283221337 * thread_id + int('p' + 'l' + 'z' + 'h' + 'e' + 'l' + 'p'));
        while (true) {
            n = prefix + cycle_nodes;
            assignment.assign(n, -1);
            edges.assign(n, {});
            vector<int> parent(prefix);
            for (int i = 1; i < prefix; i++) {
                parent[i] = rnd() % i;
                edges[parent[i]].push_back(i);
            }
            //todo add chords?
            for (int i = 0; i < cycle_len; i++) {
                assignment[prefix + i] = cyc_ass[i];
                //assignment[prefix + i] = rnd() % (players - 1);
            }
            for (auto [v, u] : cyc_edges)
                edges[prefix + v].push_back(prefix + u);
            for (int i = 0; i < prefix; i++) {
                for (int j = i + 1; j < prefix; j++)
                    if (parent[j] != i) {
                        if (rnd() % (prefix + 2) == 0)
                            edges[i].push_back(j);
                    }
                for (int j = 0; j < cycle_len; j++) {
                    if (rnd() % ((prefix + 1)) == 0)
                        edges[i].push_back(prefix + j);
                }
            }
            int additional = 0;
            for (int i = 0; i < prefix; i++) {
                if (edges[i].empty()) {
                    edges[i].push_back(i + 1 + rnd() % (prefix + cycle_len - 1 - i));
                }
                if (additional < 3 && rnd() % max(1, prefix) == 0) {
                    additional++;
                    edges.emplace_back();
                    assignment.push_back(-1);
                    edges[i].push_back(n);
                    n++;
                }
                if (!edges[i].empty()) {
                    assignment[i] = rnd() % players;
                } else {
                    assert(0);
                    assignment[i] = -1;
                }
            }
            int E_COUNT = 0;
            for (int i = 0; i < n; i++)
                E_COUNT += edges[i].size();
            if (E_COUNT > 100)
                continue;
            {
                int Z = make_checks();
                if (Z != 0) {
                    cerr << "Failed check with TOT=" + to_string(Z) + " at thread " + to_string(thread_id) + "\n";
                    continue;
                }
            }
            precalc(false);
            precheck = 0;
            if (solve(false)) {
                cerr << "Found non nash at " + to_string(IT) + "\n";
                require_exit = true;
            }
            if (IT % 50 == 0)
                cerr << "Done " + to_string(IT) + " at thread " + to_string(thread_id) + "\n";
            IT++;
            if (require_exit)
                break;
        }
        //cerr << "There is " << good_count << " non nash equilibriums in total\n";
        //cerr << "ratio is " << good_count / (double)(min(IT, r_start) - max(0, l_start)) << '\n';
    }


    pair<int, int> term_cyc_seg;

    void gen_tree_cycle() {
        ifstream cin("../a.in");
        int prefix, cycle_len, term_count;
        cin >> prefix >> cycle_len >> players >> term_count;
        if (thread_id == 0)
            cerr << "n p s " << prefix << ' ' << cycle_len << ' ' << players << '\n';
        start = 0;
        int IT = 0;
        //auto perm_n = generate_permutations(n);
        //auto perm_c = generate_permutations(cycle_len);
        mt19937 rnd(12540 * thread_id + 223);
        while (true) {
            n = prefix + cycle_len + term_count;
            term_cyc_seg = {prefix + cycle_len, n};
            assignment.assign(n, -1);
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
                assignment[prefix + i] = rnd() % players;
                if (rnd() % ((cycle_len + 1)) == 0)
                    edges[prefix + i].push_back(prefix + cycle_len + rnd() % term_count);
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
                if (i > 0 && additional < 2 && rnd() % max(1, prefix) == 0) {
                    additional++;
                    edges.emplace_back();
                    assignment.push_back(-1);
                    edges[i].push_back(n);
                    n++;
                }
                if (!edges[i].empty()) {
                    assignment[i] = rnd() % players;
                } else {
                    assert(0);
                    assignment[i] = -1;
                }
            }
            {
                int Z = make_checks();
                if (Z != 0) {
                    cerr << "Failed check with TOT=" + to_string(Z) + " at thread " + to_string(thread_id) + "\n";
                    continue;
                }
            }
            precalc(false);
            precheck = 0;
            if (solve(false)) {
                cerr << "Found non nash at " + to_string(IT) + "\n";
                require_exit = true;
            }
            if (IT % 50 == 0)
                cerr << "Done " + to_string(IT) + " at thread " + to_string(thread_id) + "\n";
            IT++;
            if (require_exit)
                break;
        }
        //cerr << "There is " << good_count << " non nash equilibriums in total\n";
        //cerr << "ratio is " << good_count / (double)(min(IT, r_start) - max(0, l_start)) << '\n';
    }

    void brute_input(int l_start, int r_start) {
        ifstream cin("../a.in");
        cin >> n >> players >> start;
        cerr << "n p s " << n << ' ' << players << ' ' << start << '\n';
        //cerr << n << ' ' << players << ' ' <<  start << '\n';
        //exit(0);
        start--;
        assignment.resize(n);
        for (int &i : assignment) {
            cin >> i;
            i--;
            assert(-2 <= i && i < players);
        }
        vector<pair<int, int>> e1, e2;
        {
            int sz;
            cin >> sz;
            e1.resize(sz);
            for (auto &[u, v] : e1) {
                cin >> u >> v;
                assert(1 <= v && v <= n);
                assert(1 <= u && u <= n);
            }
        }
        {
            int sz;
            cin >> sz;
            e2.resize(sz);
            for (auto &[u, v] : e2) {
                cin >> u >> v;
                assert(1 <= v && v <= n);
                assert(1 <= u && u <= n);
            }
        }
        cerr << "We have " << e1.size() << " Guaranteed edges and " << e2.size() << " Maybe edges\n";
        vector<vector<int>> e(n);
        for (auto [v, u] : e1)
            e[v - 1].push_back(u - 1);
        vector<int> ass = assignment;
        vector<pair<vector<int>, vector<vector<int>>>> how_much;
        int IT = 0;
        bool solve_inplace = true;
        int good_count = 0;
        int expecting;
        auto update_expecting = [&]() {
            expecting = X.fetch_add(1);
        };
        update_expecting();
        function<void(int)> brute2 = [&](int id) {
            if (id == n) {
                counter++;
                //if (counter % 10 == 0)
                //    cerr << "Done " << counter << '\n';
                edges = e;
                make_checks();
                precalc(false);
                precheck = 1;
                if (solve(false)) {
                    if (IT == expecting) {
                        if (solve_inplace) {
                            precheck = 0;
                            if (solve(false)) {
                                good_count++;
                                cerr << "Found non nash at " << IT << '\n';
                                require_exit = true;
                            }
                            if (IT % 500 == 0)
                                cerr << "Done " + to_string(IT) + "\n";
                        } else {
                            how_much.emplace_back(assignment, edges);
                        }
                        update_expecting();
                        if (IT % 1000000 == 0)
                            cerr << "Bruted " << IT << '\n';
                    }
                    IT++;
                    if (IT == r_start) {
                        cerr << "Done bruting at " << IT << '\n';
                        exit(0);
                    }
                }
                return;
            }
            assignment[id] = ass[id];
            if (assignment[id] >= 0) {
                if (e[id].empty())
                    return;
                brute2(id + 1);
            } else if (assignment[id] == -1) {
                if (!e[id].empty())
                    return;
                brute2(id + 1);
            } else {
                assert(ass[id] == -2);
                if (e[id].empty())
                    return;
                for (int i = 0; i < players; i++) {
                    assignment[id] = i;
                    brute2(id + 1);
                    if (require_exit)
                        return;
                }
            }
        };
        function<void(int, int)> brute1 = [&](int id, int able_to_add) {
            if (id == e2.size()) {
                assignment = ass;
                brute2(0);
                return;
            }
            brute1(id + 1, able_to_add);
            if (require_exit)
                return;
            if (able_to_add > 0) {
                auto[v, u] = e2[id];
                e[v - 1].push_back(u - 1);
                brute1(id + 1, able_to_add - 1);
                e[v - 1].pop_back();
            }
        };
        for (int x = 10; x <= e2.size(); x++) {
            brute1(0, x);
            if (require_exit)
                return;
            cerr << "Done searching with " << x << " maybe edges\n";
        }
        cerr << "Total of " << IT << " of something\n";
        cerr << "Done bruting all\n";

        //for (int it = l_start; it < min((int)how_much.size(), r_start); it++) {
        if (!solve_inplace) {
            mt19937 rnd(223);
            shuffle(how_much.begin(), how_much.end(), rnd);
            IT = max(0, l_start);
            for (auto[ass, ed] : how_much) {
                assignment = ass;
                edges = ed;
                precalc(false);
                precheck = 0;
                if (solve(false)) {
                    good_count++;
                    cerr << "Found non nash at " << IT << '\n';
                    //exit(0);
                }
                if (IT % 200 == 0) {
                    cerr << "Done " + to_string(IT) + "\n";
                }
                IT++;
            }
        }
        cerr << "There is " << good_count << " non nash equilibriums in total\n";
        cerr << "ratio is " << good_count / (double)(min(IT, r_start) - max(0, l_start)) << '\n';
    }

    void input(int n_, int players_, int start_, vector<int> const &assignment_, vector<vector<int>> const &edges_) {
        n = n_;
        players = players_;
        start = start_;
        assignment = assignment_;
        edges = edges_;
        make_checks();
        precalc();
    }
    void read() {
        cin >> n >> players;
        cin >> start;
        start--;
        edges.resize(n);
        assignment.resize(n);
        for (int i = 0; i < n; i++) {
            cin >> assignment[i];
            assert(assignment[i] >= 0 && assignment[i] <= players);
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
        cerr << "n players\n";
        cerr << "source\n";
        cerr << "# each a_i in [0; players], 0 means terminal\n";
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
        cerr << "4 2\n"
                "1\n"
                "1 2 2 0\n"
                "\n"
                "1\n"
                "4\n"
                "1 2\n"
                "1 3\n"
                "2 4\n"
                "3 4\n";
        cerr << "\nSame example:\n";
        cerr << "4 2\n"
                "1\n"
                "1 2 2 0\n"
                "\n"
                "2\n"
                "0 1 1 0\n"
                "0 0 0 1\n"
                "0 0 0 1\n"
                "0 0 0 0\n";

    }

    /*int encode(vector<int> const &id, vector<int> const &val) {
        assert(id.size() == val.size());
        int res = 0;
        for (int i = 0; i < val.size(); i++) {
            assert(id[i] < val[i]);
            res = res * val[i] + id[i];
        }
        return res;
    }

    vector<int> decode(int w, vector<int> const &val) {
        vector<int> res(val.size());
        for (int i = int(val.size()) - 1; i >= 0; i--) {
            res[i] = w % val[i];
            w /= val[i];
        }
        return res;
    }*/

    //encode strategy ids except for missing player
    //missing might be -1
    int encode(vector<int> const &id, int missing) {
        int res = 0;
        for (int i = 0; i < players; i++)
            if (i != missing) {
                assert(id[i] < strategy[i].size());
                res = res * strategy[i].size() + id[i];
            }
        return res;
    }

    //decode strategy ids except for missing player
    //missing player would have -1 strategy index
    //should be encoded with encode function
    vector<int> decode(int w, int missing) {
        vector<int> res(players);
        for (int i = players - 1; i >= 0; i--) {
            if (i == missing) {
                res[i] = -1;
            } else {
                res[i] = w % strategy[i].size();
                w /= strategy[i].size();
            }
        }
        return res;
    }

    //total number of strategies for all players except missing
    //missing might be -1
    int total_encode(int missing) {
        int res = 1;
        for (int i = 0; i < players; i++) {
            if (i != missing) {
                res *= strategy[i].size();
            }
        }
        return res;
    }
};
