# rhp_root_bound

Radius containing every right-half-plane root of $D(s) + N(s)e^{-sT} = 0$,
for every delay $T \ge 0$.

## Syntax

```text
R = rhp_root_bound(N, D)
```

## Description

If $\mathrm{Re}\,s \ge 0$ then $|e^{-sT}| \le 1$, so a root satisfies
$|D(s)| \le |N(s)|$. With $D$ normalized to be monic,
$|D(s)| \ge r^m - \sum|d_k|r^k$ and $|N(s)| \le \sum |n_k| r^k$ for
$r = |s|$. `R` is the positive root of
$r^m - \sum_{k<m}|d_k|r^k - \sum_k |n_k|r^k = 0$; every root with
$\mathrm{Re}\,s \ge 0$ satisfies $|s| \le R$. Requires $\deg N < \deg D$.

## Example

```matlab
R = rhp_root_bound(2, [1 1])        % 3: |s+1| <= 2 implies |s| <= 3
```

## See also

[`unstable_root_count`](unstable_root_count.md).
