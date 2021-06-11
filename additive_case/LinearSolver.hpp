#pragma once
//#pragma GCC optimize("Ofast")
//#pragma GCC optimize("fast-math")
#include "Solver.hpp"
#include <bits/stdc++.h>
using namespace std;

struct LP {
    static const double eps;

    static bool eq(double x, double y) {
        return abs(x - y) < eps;
    }

    static bool ls(double x, double y) {
        return x + eps < y;
    }

    static vector<double> simplex(vector<vector<double> > a, bool certificate) {
        int n = (int)a.size() - 1;
        int m = (int)a[0].size() - 1;
        vector<int> left(n + 1), up(m + 1);
        iota(up.begin(), up.end(), 0);
        iota(left.begin(), left.end(), m);
        auto pivot = [&](int x, int y) {
            swap(left[x], up[y]);
            double k = a[x][y];
            a[x][y] = 1;
            vector<int> vct;
            for (int j = 0; j <= m; j++) {
                a[x][j] /= k;
                if (!eq(a[x][j], 0)) vct.push_back(j);
            }
            for (int i = 0; i <= n; i++) {
                if (eq(a[i][y], 0) || i == x) continue;
                k = a[i][y];
                a[i][y] = 0;
                for (int j : vct) a[i][j] -= k * a[x][j];
            }
        };
        while (true) {
            int x = -1;
            for (int i = 1; i <= n; i++) if (ls(a[i][0], 0) && (x == -1 || a[i][0] < a[x][0])) x = i;
            if (x == -1) break;
            int y = -1;
            for (int j = 1; j <= m; j++) if (ls(a[x][j], 0) && (y == -1 || a[x][j] < a[x][y])) y = j;
            if (y == -1) return {};
            if (y == -1) assert(0); // infeasible
            pivot(x, y);
        }
        if (!certificate)
            return {1};
        while (true) {
            int y = -1;
            for (int j = 1; j <= m; j++) if (ls(0, a[0][j]) && (y == -1 || a[0][j] > a[0][y])) y = j;
            if (y == -1) break;
            int x = -1;
            for (int i = 1; i <= n; i++) if (ls(0, a[i][y]) && (x == -1 || a[i][0] / a[i][y] < a[x][0] / a[x][y])) x = i;
            if (x == -1) assert(0); // unbounded
            pivot(x, y);
        }
        vector<double> ans(m + 1);
        for (int i = 1; i <= n; i++) if (left[i] <= m) ans[left[i]] = a[i][0];
        ans[0] = -a[0][0];
        return ans;
    }
    // j=1..m: x[j]>=0
    // i=1..n: sum(j=1..m) A[i][j]*x[j] <= A[i][0]
    // max sum(j=1..m) A[0][j]*x[j]
    // res[0] is answer
    // res[1..m] is certificate

    vector<pair<vector<int>, int>> q;

    vector<int> added;

    void add(vector<pair<vector<int>, int>> const &s) {
        for (auto const &i : s)
            q.push_back(i);
        added.push_back(s.size());
    }

    void remove() {
        for (int i = 0; i < added.back(); i++)
            q.pop_back();

        added.pop_back();
    }
//#define NEGATIVE

#ifdef NEGATIVE
    vector<double> get(bool certificate) {
        vector<vector<double>> r;
        if (q.empty()) {
            assert(!certificate);
            return {1};
        }
        int n = q[0].first.size();
        r.emplace_back(2 * n + 1);
        for (int i = 0; i < n; i++) {
            r[0][i + 1] = -1;
            r[0][n + i + 1] = 1;
        }
        r[0][0] = 0;
        if (certificate) {
            r.push_back(r[0]);
            r.back()[0] = 10;
        }
        for (auto const &i : q) {
            r.emplace_back(2 * n + 1);
            r.back()[0] = -i.second;
            for (int j = 0; j < n; j++) {
                r.back()[j + 1] = -i.first[j];
                r.back()[j + 1 + n] = i.first[j];
            }
        }
        return simplex(r, certificate);
    }
#else

    vector<double> get(bool certificate) {
        vector<vector<double>> r;
        if (q.empty()) {
            assert(!certificate);
            return {1};
        }
        r.emplace_back(q[0].first.size() + 1, -1);
        r[0][0] = 0;
        for (auto const &i : q) {
            r.emplace_back();
            r.back().push_back(-i.second);
            for (int j : i.first)
                r.back().push_back(-j);
        }
        return simplex(r, certificate);
    }
#endif
};

const double LP::eps = 1e-9;

vector<int> operator-(vector<int> const &a, vector<int> const &b) {
    assert(a.size() == b.size());
    vector<int> c(a);
    for (int i = 0; i < (int)a.size(); i++)
        c[i] -= b[i];
    return c;
}

mt19937 reee(231323);
// For each strategy of player x, chooses best strategy for player y,
// then looks for solution via linear solver
struct LinearSolver : Solver {

    LinearSolver() = default;

    vector<vector<int>> number;
    vector<vector<vector<int>>> paths;

    vector<vector<int>> first[2];
    vector<vector<vector<int>>> first_paths[2];

    int N, M;


    vector<int> vis;

    LP solvers[2];

    int step_size = 3;

    int biggest = 0;

    int total = 0;

    int random_seed;

    template<typename T, typename U>
    void shuffle(vector<T> &a, vector<U> &b) {
        //return;
        assert(a.size() == b.size());
        mt19937 rnd(random_seed);
        for (int i = 0; i < a.size(); i++) {
            int j = rnd() % (i + 1);
            swap(a[i], a[j]);
            swap(b[i], b[j]);
        }
    }
int cut = 0;
    bool go(int it = 0, int tp = 0) {
        if (biggest < it || (biggest == it && (reee() & 31) == 0)) {
            cut++;
            //cerr << it << ' ' << cut << '\n';
        }
        biggest = max(biggest, it);
        //if (it == biggest)
        //    cerr << it << '\n';
        bool checked = false;
        auto check = [&]() {
            //return true;
            if (checked)
                return true;
            return checked = (!solvers[0].get(false).empty() && !solvers[1].get(false).empty());
        };

        if ((step_size != -1 && tp == 0 && it % step_size == 0) || it == max(N, M))
            if (!check())
                return false;

        if (it == max(N, M)) {
            //total++;
            //cerr << total << '\n';
            //return true;
            cerr << "FOUND\n";
            for (int i = 0; i < n; i++)
                cerr << assignment[i] << ' ';
            cerr << '\n';
            vector<double> res[2];
            for (int t = 0; t < 2; t++) {
                res[t] = solvers[t].get(true);
                assert(!res[t].empty());
            }
            print_time();
            cout << fixed;
            cout.precision(5);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < (int)edges[i].size(); j++) {
                    int id = edges_id[i][j];
                    cerr << i + 1 << ' ' << edges[i][j] + 1 << ' ' <<
                                                                   #ifndef NEGATIVE
                    res[1][id + 1] << ' ' << res[0][id + 1] << '\n';
                                                                   #else
                    res[1][id + 1] - res[1][m + id + 1] << ' ' << res[0][id + 1] - res[0][m + id + 1] << '\n';
                                                                   #endif
                }
            }
            exit(0);
        }
        if (it == 8) {
            int kek = 1;
        }
        if (it >= (int)strategy[tp].size() || first[tp][it].empty())
            return go(it + tp, !tp);

        int req_id = -1;

        /*
        for (int i = 0; i < (int)first[tp][it].size(); i++)
            if (vis[first[tp][it][i]] == tp) {
                if (req_id == -1)
                    req_id = i;
                else
                    return false;
            }
        */

        int curs = 0;
        auto try_do = [&](int i) {
            int id = first[tp][it][i];
            if (vis[id] == -1 || vis[id] == tp) {
                int old = vis[id];
                vis[id] = tp;

                vector<pair<vector<int>, int>> now;
                auto cur = first_paths[tp][it][i];
                for (int j = 0; j < (int)first[tp][it].size(); j++)
                    if (i != j)
                    //if (vis[first[tp][it][j]] == !tp)
                        now.emplace_back(first_paths[tp][it][j] - cur, 1);
                    //asd
                    //else if (vis[first[tp][it][j]] == -1)
                    //    now.emplace_back(first_paths[tp][it][j] - cur, 0);

                solvers[tp].add(now);
                bool res = go(it + tp, !tp);
                vis[id] = old;
                solvers[tp].remove();
                if (!res) {
                    curs++;
                } else {
                    curs = 0;
                }
                //5 works
                if (curs == 5)
                    return false;
                if (!res && !check())
                    return false;
            }
            return true;
        };

        if (req_id != -1) {
            if (!try_do(req_id))
                return false;
        } else {
            //shuffle(first[tp][it], first_paths[tp][it]);
            for (int i = 0; i < (int)first[tp][it].size(); i++)
                if (!try_do(i))
                    return false;
        }

        return true;
    }


    vector<int> get_score(int i1, int j1) {
        //ClpSimplex model;
        //model.primal(1);
        vector<int> out(n, -1);
        for (int t = 0; t < 2; t++) {
            int id = vector<int>{i1, j1}[t];
            for (int i = 0; i < (int)assigned[t].size(); i++)
                out[assigned[t][i]] = strategy_id[t][id][i];
        }
        vector<int> res(m);
        vector<bool> visited(n);
        int v = start;
        while (!visited[v] && v != finish) {
            visited[v] = true;
            int id = out[v];
            assert(id >= 0);
            res[edges_id[v][id]] = 1;
            v = edges[v][id];
        }
        if (v != finish) {
            return {};
        }
        return res;
    }


    void precalc() {
        {

            step_size = 3;

            biggest = 0;

            total = 0;

            cut = 0;

            number.clear();
            paths.clear();
            vis.clear();
            for (int t = 0; t < 2; t++) {
                first[t].clear();
                first_paths[t].clear();
                solvers[t].q.clear();
                solvers[t].added.clear();
            }
        }
        /*for (int t = 0; t < 2; t++) {
            cout << "player " << t + 1 << '\n';
            for (int i = 0; i < strategy[t].size(); i++) {
                cout << i << '\n';
                for (int j = 0; j < assigned[t].size(); j++)
                    cout << assigned[t][j] << ' ' << strategy[t][i][j] << ' ' << strategy_id[t][i][j] << '\n';
            }
        }
        exit(0);*/
        //for (int i = 0; i < n; i++)
        //    for (int j = 0; j < (int)edges[i].size(); j++)
        //q        cerr << i + 1 << ' ' << edges[i][j] + 1 << ' ' << edges_id[i][j] << '\n';

        if (random_seed == -1)
            random_seed = chrono::high_resolution_clock::now().time_since_epoch().count();

        if (false) {
            for (int t = 0; t < 2; t++)
                shuffle(first[t], first_paths[t]);
        }

        vector<long long> random_values(m);
        {
            mt19937_64 rnd(random_seed);
            for (int i = 0; i < m; i++)
                random_values[i] = rnd();
        }

        //strategy[1].pop_back();
        //strategy_id[1].pop_back();

        N = strategy[0].size();
        M = strategy[1].size();

        number.resize(N, vector<int>(M, -1));
        paths.resize(N, vector<vector<int>>(M));
        map<long long, int> mp;
        mp[0] = -1;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++) {
                auto g = get_score(i, j);
                paths[i][j] = g;
                long long cur = 0;
                for (int x = 0; x < (int)g.size(); x++)
                    if (g[x] == 1)
                        cur += random_values[x];
                if (mp.count(cur) == 0) {
                    int id = mp.size();
                    mp[cur] = id;
                }
                number[i][j] = mp[cur];
            }
        cerr << "Total of " << mp.size() << " paths\n";
        for (int t = 0; t < 2; t++) {
            int sz = strategy[t].size();
            first[t].resize(sz);
            first_paths[t].resize(sz);
            for (int i = 0; i < sz; i++) {
                set<int> exist;
                for (int j = 0; j < (int)strategy[!t].size(); j++) {
                    int x = i, y = j;
                    if (t)
                        swap(x, y);
                    int num = number[x][y];
                    if (num != -1 && exist.find(num) == exist.end()) {
                        first[t][i].push_back(num);
                        first_paths[t][i].push_back(paths[x][y]);
                        exist.insert(num);
                    }
                }
            }
        }
        vis.resize(mp.size(), -1);
    }

    clock_t start_time;

    void print_time() {
        cerr << "Consumed " << (clock() - start_time) / (double)CLOCKS_PER_SEC << '\n';
    }

    void solve() {
        random_seed = 51231;
        start_time = clock();
        precalc();

        go();
        //cerr << "Total of " << total << '\n';

        ..cerr << "Not found\n";
        //print_time();
    }
};
