# Lattice Boltzmann Method D2Q9
Implementation in Fortran by Bartosz Mierzwa prepared for *High-Performance Computing: Fluid Mechanics with Fortran* 
course at Albrecht-Ludwig's University Freiburg based on lecture by Andreas Greiner.

This repository consist of number of branches:
* **master**: contains the most universal LBM d2q9 solver
  * it is highly configurable, allowing for simulating many types of fluids movement. It offers e.g.:
    * 1 external force: pressure
    * 3 initial field conditions: constant, sine wave, gaussian (applicable to both velocity and density fields)
    * 2 types of boundry conditions: periodic and bounceback
    * optional simulation of moving wall
    * boundry conditions setup on every wall
    * grid size, time steps and relaxation parameter $\omega$ configuration
  * code is parallelized using coarray feature
  * program allows for results output logging. It offers:
    * 4 loggers: field, profile_x, profile_y and point
    * number of frames configuration
    * debug info on loger jobs
  * 7 experiments families are saved in config directory:
    * gaussian-wave-decay - simulates gaussian density pulse in a box and a tube
    * density-wave-decay - simulates density sine wave decay
    * shear-wave decay - simulates velocity shear wave decay + viscosity measures
    * couette-flow - simulates couette flow in a tube
    * poiseuille-flow - simulates poiseuille flow in a tube with pressure gradient
    * sliding-lid - simulates lid driven cavity flow in a box
    * moving-walls - more simulations with moving walls
  * results presentations based on experiment made on BW UNI CLUSTER
    * results are shown in eval/results.ipynb notebook
    * experiment data are not available due to repository size restrictions
* **optim** - branch family where few optimization were tested for optimal sliding-lid simulation:
  * there are 4 sub-branches:
    * minimal - limits code to minimum necessary to run the experiment and fuses code collision into single subroutine
    * pull-streaming - replaces push with pull streaming strategy (sequential writes instead reads during streaming)
    * fusion - fusion of streaming and collision subroutines
    * spaghetti - attempts to limit overhead by collapsing code into spaghetti
  * notes:
    * parameters must be hard-coded before execution
    * each experiment outputs: 
      * last lbm lattice state
      * statistics: num of tasks, image number, MLUPS, total execution time, collision time, streaming time, collision time
      * statistics do not take into account time spent on IO operations
* **scaling** - branch family where few scaling experiments were tested:
  * there are 4 sub-branches:
    * Re1000-L300 - evaluates scaling on ice lakes cpus with Re number 1000 on grid size 300 x 300
    * Re1000-L1536 - evaluates scaling on ice lakes cpus with Re number 1000 on grid size 1536 x 1536
    * Re1000-L1536-spread - evaluates scaling on ice lakes cpus with Re number 1000 on grid size 1536 x 1536 where each node usage is limited to 32 tasks
    * Re1000-L1536-cpu - evaluates scaling on BW UNI CLUSTER 96-task's cpus with Re number 1000 on grid size 1536 x 1536
  * notes: 
    * evaluation graphs are presented in eval\evaluation.ipynb
    * bear in mind that on the heatmaps the Y direction is reversed, so top is SOUTH and bottom is NORTH
  