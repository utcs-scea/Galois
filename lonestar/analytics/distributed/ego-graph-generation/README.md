
```
./ego-graph-generation-dist -startNode=<int Global ID, default to 0> <inputGraph.gr>
```

# Algorithm

## All-source ego-graph generation (ASEG)

```
L: an array of integers specifying number of nodes per level

NodeData:
- NL: list of sets, each set containing global node ids for each level
- EL: list of sets, same as NL but for edges

// Initialization
for each node on each host:
    NL <- [{}, ..., {}]
    EL <- [{}, ..., {}]

// First iteration for level 0
for each node src on each host:
    for each (outgoing) neighbor dst of src:
        if size of NL[0] reaches L[0]:  // check first in case level limit is 0
            break
        if dst == src:  // ignore self loops
            continue
        add global id of dst to NL[0] of src
        add (dst, dst), (src, dst), (dst, src) to EL[0] of src
sync NL and EL among masters and mirrors

// Iterations for remaining levels
for level l in 1:size(L):
    for each node src on each host:
        for each node dst0 in NL[0] of src:
            for each node dst in NL[l - 1] of dst0:
                // Rational: dst0 is at level 0 for src, dst is at level (l - 1)
                //  for dst0, so dst0 can be at level l for src
                // Note: dst0 is an immediate neighbor for src and thus is
                //  accessible from src, but dst can be a node on another host
                //  and thus may not be accessible from either src or dst0 (when
                //  dst0 is a mirror node under oec policy)
                if size of NL[l] reaches L[l]:
                    break
                if dst == src or dst in NL of src at any level:
                    continue
                add dst to NL[l] of src
            for each edge e in EL[l - 1]:
                if either end of e is in NL of src at any level:
                    add e to EL[l] of src
    sync NL and EL among masters and mirrors
```

### Issue

Suppose we have a tree-like graph
```
         / -> 3
    -> 1
  /      \ -> 4
0
  \      / -> 5
    -> 2
         \ -> 6
```
and we partition it to 2 hosts with oec
```
Host 1          Host2
              |        / -> 3(2)
    -> 1(1,m)   1(0,M)
  /           |        \ -> 4(3)
0(0)
  \           |        / -> 5(4)
    -> 2(2,m)   2(1,M)
                       \ -> 6(5)
```
where local ids are in parenthesis; M is master and m is mirror.

Assume each node has its children as level 0 nodes. When choosing nodes for level 1 for node 0 on host 1,
let's assume we choose global 5, 6, 4 on host 2. These global node ids are gathered from immediate neighbors 1 and 2.
However, if there is an edge between 4 and 5 (which is 3 and 4 on host 2), node 0 can hardly capture this edge in its ego-graph
because vertex programming model only considers local/immediate neighborhood.