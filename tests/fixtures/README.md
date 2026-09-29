`kernel_v080.mat` was generated with the unmodified 0.8.0 code at 9b70d02,
in MATLAB R2026a, before the shared-grid implementation:

```matlab
t=(0:.1:1).';
scalar=ckernel(ndpair(1,[1 1]),t);
matrix=ckernel(cmimo({ndpair(1,[1 1]),2;0,ndpair(2,[1 2])}),t);
save('kernel_v080.mat','scalar','matrix','t');
```

These immutable numeric schema-1 snapshots exercise loading actual old
objects, rather than saving and loading only the current implementation.
