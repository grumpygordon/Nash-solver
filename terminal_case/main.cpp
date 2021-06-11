#undef NDEBUG
//#pragma GCC optimize("Ofast")

#include <bits/stdc++.h>
using namespace std;
#include "RandomSolver.hpp"
#include "LinearSolver.hpp"
#include <unistd.h>

int main(int argc, char** argv) {
#ifdef ONPC
    freopen("../a.in", "r", stdin);
    freopen("../a.out", "w", stdout);
#endif
    cin.tie(0);
    cout.tie(0);
    //RandomSolver q;
    //signal (SIGINT,my_handler);
    atomic<int> X{0};
    vector<thread> threads;
    atomic_bool require_exit{false};
    if (false) {
        freopen("../a.in", "r", stdin);
        RandomSolver q(ref(X), ref(require_exit), 0);
        q.read();
        q.solve();
        return 0;
    }
    LinearSolver q(ref(X), ref(require_exit), 0);
    //q.help();
    //q.read();
    int l = 0, r = 1e9;
    //int l = atoi(argv[1]), r = atoi(argv[2]);
    //cerr << "Will search in [" << l << "; " << r << ")\n";
    for (int it = 0; it < 7; it++) {
        threads.emplace_back([it, &X, &require_exit]() {
            LinearSolver q(ref(X), ref(require_exit), it);
            //q.help();
            //q.read();
            int l = 0, r = 1e9;
            //int l = atoi(argv[1]), r = atoi(argv[2]);
            //cerr << "Will search in [" << l << "; " << r << ")\n";
            //q.gen_tree_cycle();
            q.gen_tree();
            //q.brute_input(l, r);
            //q.solve();
        });
    }
    for (auto &i : threads)
        i.join();
    cerr << "Consumed " << (clock() / (double)CLOCKS_PER_SEC) << '\n';
}
