---
title: CF1073G Yet Another LCP Problem
date: 2024-01-21 18:23:58
tags: 
    - String
    - SAM
    - 虚树
    - DP
categories:
    - 题解
    - Codeforces
mathjax: true
---

还算有点意思。

把原字符串倒过来建一个 SAM，则我们熟知两个后缀的 LCP 就是其对应节点在 link 树上 LCA 的 `len` 值。

然后考虑 DP。这类问题的普遍套路是钦定某一个点是 LCA 然后对每个点计算贡献。我们设 $f_i$ 表示 $i$ 子树内属于集合 $A$ 的点的个数，$g_i$ 表示 $j$ 子树内属于集合 $B$ 的点的个数，显然一对合法的点的 LCA 是 $i$ 当且仅当它们在 $i$ 的不同子树中，因此，我们便有如下方程：

$$
ans = \sum_{i} (g_u [i \in A] + \sum_{j \in \operatorname{subtree}(i)} f_j \cdot (g_i - g_j)) \cdot len_i
$$

直接做会退化成 $\mathcal{O}(qn)$，但考虑到 $\sum |A|$ 及 $\sum |B|$ 的量级都为 $2 \times 10^5$，因此我们可以考虑建出虚树然后再 DP，这样复杂度就对了。

注意清空虚树的邻接表以及 $f$ 与 $g$。

```cpp
#include <bits/stdc++.h>

using i64 = long long;
using u32 = unsigned;

template <typename T = int>
inline T read() {
  T res;
  std::cin >> res;
  return res;
}

constexpr int N = 4e5 + 50;

int n, m;
char s[N];

int pos[N];

struct SAM {
  int ch[N][26], link[N], len[N], cnt = 1, lst = 1;
  int expand(int d) {
    int cur = ++cnt;
    len[cur] = len[lst] + 1;
    int p = lst;
    for (; p && !ch[p][d]; p = link[p]) ch[p][d] = cur;
    if (!p) {
      link[cur] = 1;
    } else {
      int q = ch[p][d];
      if (len[p] + 1 == len[q]) {
        link[cur] = q;
      } else {
        int tmp = ++cnt;
        len[tmp] = len[p] + 1, link[tmp] = link[q];
        memcpy(ch[tmp], ch[q], sizeof(ch[q]));
        for (; p && ch[p][d] == q; p = link[p]) ch[p][d] = tmp;
        link[cur] = link[q] = tmp;
      }
    }
    return lst = cur;
  }
} sam;

std::vector<int> adj[N];

namespace LCA {

int dfn[N], clk;
int st[20][N];

void dfs(int u, int frm) {
  st[0][dfn[u] = ++clk] = frm;
  for (auto v : adj[u]) {
    if (v == frm) continue;
    dfs(v, u);
  }  
}

inline int cmp(int u, int v) { return dfn[u] < dfn[v] ? u : v; }

void init() {
  dfs(1, 0);
  for (int i = 1; i <= std::__lg(clk); i++) {
    for (int j = 1; j <= clk - (1 << i) + 1; j++) {
      st[i][j] = cmp(st[i - 1][j], st[i - 1][j + (1 << (i - 1))]);
    }
  }
}

inline int queryLCA(int u, int v) {
  if (u == v) return u;
  u = dfn[u], v = dfn[v];
  if (u > v) std::swap(u, v); 
  int d = std::__lg(v - u++);
  return cmp(st[d][u], st[d][v - (1 << d) + 1]);
}

}  // namespace LCA

using LCA::queryLCA, LCA::dfn;

std::vector<int> vadj[N];
int stk[N], top;
i64 f[N], g[N], ans;

void dp(int u) {
  int hav = f[u], len = sam.len[u];  // hav == 1 means there is a point in u 
  for (auto v : vadj[u]) dp(v), f[u] += f[v], g[u] += g[v];
  for (auto v : vadj[u]) ans += f[v] * (g[u] - g[v]) * len;
  ans += hav * g[u] * len;
}

void solve() {
  std::cin >> n >> m >> (s + 1);
  for (int i = n; i >= 1; i--) pos[i] = sam.expand(s[i] - 'a');

  for (int i = 2; i <= sam.cnt; i++) adj[sam.link[i]].push_back(i);
  LCA::init();

  for (int k, l, u; m; m--) {
    std::cin >> k >> l;
    
    std::vector<int> t;
    while (k--) t.push_back(pos[u = read()]), f[pos[u]]++;
    while (l--) t.push_back(pos[u = read()]), g[pos[u]]++;

    std::sort(t.begin(), t.end(), [&](int lhs, int rhs) {
      return dfn[lhs] < dfn[rhs];
    });
    t.erase(std::unique(t.begin(), t.end()), t.end());
    u32 siz = t.size();
    for (u32 i = 1; i < siz; i++) t.push_back(queryLCA(t[i - 1], t[i]));
    std::sort(t.begin(), t.end(), [&](int lhs, int rhs) {
      return dfn[lhs] < dfn[rhs];
    });
    t.erase(std::unique(t.begin(), t.end()), t.end());

    stk[top = 1] = t[0];
    for (u32 i = 1; i < t.size(); i++) {
      while (queryLCA(t[i], stk[top]) != stk[top]) top--;
      vadj[stk[top]].push_back(t[i]), stk[++top] = t[i];
    }

    ans = 0, dp(t[0]);
    std::cout << ans << "\n";

    for (auto i : t) vadj[i].clear(), f[i] = g[i] = 0;
  }
}

int main() {
  std::cin.tie(nullptr)->sync_with_stdio(false);

  int t = 1;
  // std::cin >> t;
  while (t--) solve();

  return 0;
}
```
