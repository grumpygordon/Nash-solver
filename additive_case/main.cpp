#include <bits/stdc++.h>

using namespace std;
//#include "RandomSolver.hpp"
#include "LinearSolver.hpp"

int main() {
    freopen("../a.in", "r", stdin);
    //freopen("../a.out", "w", stdout);
    cin.tie(0);
    cout.tie(0);
    //ClpSimplex s;
    LinearSolver q;
    //q.help();
    q.read();
    q.solve();
    //q.gen_tree();
    //q.gen_tree_cycle();
}
