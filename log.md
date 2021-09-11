# 2021.9.1

`ising-chain-prototype\2021-9-1-old` is created using previous codes.

By running `transverse_ising-1d-with-ani-ising.jl` two images about magnetization are created.

`ising-chain-prototype\2021-9-1-old` is then rewritten into a more modularized version, stored in `ising-chain-prototype\2021-9-1-revised`.

`ising-chain-prototype\2021-9-1-revised` is later used to created `ising-field-coupling-prototype\2021-9-1`.

# 2021.9.7

`ising-chain-prototype\2021-9-1-revised` is copied to `src\pure-z2`, as the prototype of the final version
of pure $\mathbb{Z}_2$ gauge field simulation.

`src\pure-z2\magnetization.PNG` and `src\pure-z2\magnetic-susceptibility.PNG` are created using `magnetic-ordering.jl`.
Console output:

> Progress: 100%|██████████████████████████████████████████████████████████████████████| Time: 0:20:04    

`src\pure-z2\wilson-loop.jl` is created. A rough test indicates that it works well, since systems predicted to be
deconfined obeys the perimeter law, while systems predicted to be confined obeys the area law.

Further tests may be needed, to find out:
- why the transition point seem to be $h = 1.2$ instead of $h = 1$
- what is the best choice of $T$
- how to simulate the finite temperature case

# 2021.9.8 

Add introduction of orthogonal metals into `docs/note.tex`.

`phase-diagram-transverse-ising-metropolis.PNG` and `phase-diagram-transverse-ising-wolff.PNG` are copied from `visconfs` to `docs/phase`.

Deprecated `fermion-coupling-prototype/hopping-hamiltonian.jl`. The file is moved into `deprecated`.

# 2021.9.10

Rename `src/pure-z2` to `src/pure-z2-mapped-to-spin-chain`, because it seems that with a strong transverse field, we cannot simply map a $\mathbb{Z}_2$ gauge field into a bundle of spin chains.

Create `analytical` folder for calculating and plotting analytical results.

# 2021.9.11

`2d-tfim-to-3d-cim.jl` and `2d-tfim-to-3d-cim-prototype.jl` are moved into `analytical/jxy-jtau-ratio`.

`2d-tfim-to-3d-cim.jl` is used to create `ratio-h=1.PNG` and `ratio-h=4.PNG`.
They are copied to `docs/montecarlo`. 