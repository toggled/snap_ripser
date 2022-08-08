Snap-ripser is an adaptation of Snap ([Graph library](https://github.com/snap-stanford/snap)) and ripser ([C++ utility for Rips persistence](https://github.com/Ripser/ripser)) to compute $\epsilon$-net induced lazy-witness persistence of weighted and unweighted graphs. The file `examples/effectiveness/lazywitnessph_effect.cpp` is the starting point. 

<b> Core-features </b>
- Landmark selection heuristics for constructing LW filtration on graphs (function `select_landmarks()`):
  - Greedy $\epsilon$-net (criteria=`"eps_baseline2"`)
  - Iterative $\epsilon$-net (criteria = `"epsnetring"`)
  - SPTpruning $\epsilon$-Net (criteria=`"eps_filtering"`)
  - Others:
    - node-degree based landmark selection
    - random landmark selection
    - eigen-vector centrality based landmark selection
    - maxmin landmark selection (distance-measure: weighted geodesic distance)
    - $\epsilon$-net maxmin
- Sparse distance matrix construction from Landmark-set $L$ to node-set $V$ (function `construct_landmark_distmat()`)
  
<h2> System Requirements & Compilation </h2>  

<b> Operating system</b>
- Ubuntu 18.04 LTS
  
<b>C/C++ Compiler </b>
- g++-5 (`sudo apt-get install gcc-5 g++-5`)
  - If `unable to locate package gcc-5`, run 
    - `sudo apt update` 
    - `sudo apt-get install gcc-5 g++-5`
- Symlink gcc-5 and g++-5 as gcc and g++ respectively.
  - `sudo ln -s /usr/bin/gcc-5 /usr/bin/gcc`
  - `sudo ln -s /usr/bin/g++-5 /usr/bin/g++`
  
<b>Dependencies</b>
- make (`sudo apt install make`)
- boost 1.79.0
  - Download boost: https://www.boost.org/doc/libs/1_79_0/more/getting_started/unix-variants.html 
  - Extract boost_1_79_0.tar.gz
  - `cd boost_1.79_0`
  - `./bootstrap.sh`
  - `sudo ./b2 install`

<b> Compiling snap_ripser </b>
- `cd snap_ripser`
- `make all`

<h2> Example Run </h2>

- Weighted graph: To run Greedy $\epsilon$-net with $\epsilon = 0.5$ for selecting landmarks on Celegans, and finally computing intervals up-to dimension 2.
  - `cd examples/effectiveness`
  - `./lazywitnessph_effect --format wgraph ../celegans_final.edges --heuristic eps_baseline2 --dim 2 --epsilon 0.1 --iterations 1`

- Unweighted graph:
  - `cd examples/effectiveness`
  - `./lazywitnessph_effect --format graph ../power.edge --heuristic eps_baseline2 --dim 2 --epsilon 6 --iterations 1`
  
<b> Script for multiple runs </b>

Script automate.sh varies $\epsilon$, computes landmarks, persistence intervals and finally generates 1) various statistics (e.g. Landmark selection time, Total computation time, #landmarks) in .time files and 2) dim-0, dim-1,.. barcodes in .csv files. Find these files under respective [algorithm_name]\_[dataset_name]\_[ $\epsilon$ ] folders.
- `cd examples/effectiveness`
- `bash automate.sh`

<h2> Citation </h2>

- Arafat, Naheed Anjum, Debabrota Basu, and Stéphane Bressan. "Topological Data Analysis with $\epsilon$-net Induced Lazy Witness Complex." International Conference on Database and Expert Systems Applications. Springer, Cham, 2019.
  
- Arafat, Naheed Anjum, Debabrota Basu, and Stéphane Bressan. " $\epsilon$-net Induced Lazy Witness Complexes on Graphs." arXiv preprint arXiv:2009.13071 (2020).



