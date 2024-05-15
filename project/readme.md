<div style="margin: 0 100px 0 100px;  font-family: 'Times New Roman';"> 
<h1 style="text-align: center; font-weight: normal;">TDDE65 Miniproject Report</h1>
<p style="text-align: center;">Albin Arvidsson albar556 <br>Martin Kaller marka727<p>
<p style="text-align: center;">May 15, 2024<p>

<div style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
<img id="figure1" src="./image.png"></img>
Figure 1: Particle simulation.
</div>

## 1&emsp;Introduction

We made a particule simulator with MPI (Message Passing Interface) <a style="color:inherit;" href="#reference1">[1]</a>, to verify the gas law $pV = nRT$. All particles had a radius of 1 and all collisions were perfectly elastic. No friction or other physical forces acted upon the particles. The particles moved around inside a rectangular space for simpler calcualtions, see <a style="color:inherit;" href="#figure1">figure 1</a>. MPI was used to speed up the calculations via parallelization. 

#### 1.1&emsp;Limitations
- A particle can only collide with exactly one other particle.
- A particle can not both collide and bounce on the wall
    - Can causes a particle to be outside of box wall for one timestamp

## 2&emsp;Method
We divided the box into a grid, where each process gets one cell. We use a 2 dimensional MPI layout for assigning the cell.
<div style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
<img id="figure2" src="./mpi.png"></img>
Figure 2: MPI processes in a 2D grid, with particles in each process cell. Communications of particles leaving cell visualised.
</div>
<br>


```cpp
MPI_Comm GRID_COMM_MPI;
MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &GRID_COMM_MPI);
int my_coords[2];
int my_rank;
MPI_Cart_get(GRID_COMM_MPI, 2, dims, periods, my_coords);
MPI_Comm_rank(GRID_COMM_MPI, &my_rank);
```
Each process generates `TOTAL_PARTICLES/processes` amount of random particles. Particles are stored in `std::vector` (faster than `std::forward_list`, see discussion).
At the end of one timestamp iteration, we check which particles are outside of the process' cell. If it has a neighbor who could recive them, we remove the particles from our particle list and send them to our neighbour.
```cpp
unsigned send_up_size = 0, send_left_size = 0, send_right_size = 0, send_down_size = 0;

for (auto p = particles.begin(); p != particles.end();)
{
    if (p->coords.x < box_border.x0 && left_rank != -1)
    {
        particles_send_left[send_left_size++] = p->coords;
        p = particles.erase(p);
    }
    else if (p->coords.x > box_border.x1 && right_rank != -1)
    {
       // ...
    }
    else if (p->coords.y < box_border.y0 && down_rank != -1)
    {
        // ...
    }
    else if (p->coords.y > box_border.y1 && up_rank != -1)
    {
        // ...
    }
    else
    {
        ++p;
    }
}
```
All communications are done asynchronously in MPI with `MPI_Isend` and `MPI_Irecv`
```cpp
if (up_rank != -1)
{
    // begin send and recive
    MPI_Isend(particles_send_up, send_up_size, PCORD_MPI, up_rank, 0, GRID_COMM_MPI, &req_up);
    MPI_Irecv(particles_recv_up, INIT_NO_PARTICLES, PCORD_MPI, up_rank, 0, GRID_COMM_MPI, &req_recv_up);
}
// ... left, down, right

if (up_rank != -1) {
    // wait for response to be done
    MPI_Status status;
    MPI_Wait(&req_recv_up, &status);
    int count;
    MPI_Get_count(&status, PCORD_MPI, &count);
    for (size_t i = 0; i < count; i++)
    {
        particles.push_back({recv_buf[i], false});
    }
}
// ... left, down, right

MPI_Status tmp_status;

if (up_rank != -1)
{
    // finally wait for send request to be done
    MPI_Wait(&req_up, &tmp_status);
}
// ... left, down, right

```


## 3&emsp;Debugging with DDT
<div style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
<img id="figure3" src="./ddt.png"></img>
Figure 3: DDT Debugger used to inspect sent data.
</div>
<br>

We used DDT for debugging, and gathering information during our miniproject. We encountered several issues during the implementation of the miniproject and used DDT to resolve them. For example issues with particles being sent to non-existent neighbors and just disappering. We noticed our pressure results strangely diminishing with higher timestamps. With the help of DDT it was easy to identify the issue and determine what caused it. Another use case we used DDT for was measuring the amount of sent data and comparing it between each process each step. This was made easy with DDT, see <a style="color:inherit;" href="#figure3">figure 3</a>.

## 4&emsp;Performance analysis with ITAC
<div style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
<img id="figure4" src="./itac.png"></img>
Figure 4: ITAC showing our program's executing trace.
</div>
<br>

We used ITAC during the miniproject to determine where the bottleneck in the system was. As can be seen in <a style="color:inherit;" href="#figure4">figure 4</a> we quickly determined that the computations were the heavy calculation, and that's where the optimization focus should lay for improved execution times. This differed from lab1b where the MPI communications were significantly more expensive than the actual calculations, so in lab1b we had to focus more on optimizing the MPI communications for better time performance. In the miniproject however we switched our focus to optimizing our application code, by for example improving cache locality by testing different data structures.

## 5&emsp;Results

#### 5.1&emsp;Ideal gas law

<div style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
    <img id="figure5" src="./size.svg" style="width: 100%;"></img>
    Figure 5: Relation between box size and temperature, with 32000 particles.
</div>
<br>
<div style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
    <img id="figure6" src="./particles.svg" style="width: 100%;"></img>
    Figure 6: Relation between particle amount and temperature, with size 1000.
</div>
<br>

The gas law $pV = nRT$ was verified by calculating the temperature $T$ based on the known values we got from running our program. If the program follows the gas law, then the temperature should remain constant for any combination of box size and particle amount. Our result can be seen in <a style="color:inherit;" href="#figure5">figure 5</a> and <a style="color:inherit;" href="#figure6">figure 6</a>. The program seems to follow the law quite well when the ratio of particles to size is small.

#### 5.2&emsp;Speedup

<br>
<div style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
    <img id="figure7" src="./time.svg" style="width: 100%;"></img>
    Figure 7: Relation between execution time and process amount, with 32000 particles.
</div>
<br>

<div id="table1" style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
Table 1: Relation between execution time and process amount, with 32000 particles.

<div>

| Processes | 1      | 2      | 4      | 8      | 16     | 32     | 64     |
|-----------|--------|--------|--------|--------|--------|--------|--------|
| Time (s)  | 1636.98 | 422.378 | 110.662 | 28.0934 | 8.29144 | 2.23403 | 0.595723 |

</div>

</div>

The speedup of the program from using multiple processes can be seen in <a style="color:inherit;" href="#figure5">figure 7</a> and <a style="color:inherit;" href="#table1">table 1</a>. The speedup is superlienear, with a speedup of almost fourfold when the amount of processes is doubled. 

#### 5.3&emsp;Linear vs Linked List

<br>
<div id="table2" style="display: flex; justify-content: center; align-items:center; flex-direction: column;">
Table 2: Comparison of performance of std::forward_list and std::vector, with 32000 particles.

<div>

| Processes | 8      | 16     | 32     | 64     |
|-----------|--------|--------|--------|--------|
| Time (s) std::forward_list  | 47.7115 | 12.7871 | 3.8445 | 0.998647 |
| Time (s) std::vector | 28.0934 | 8.29144 | 2.23403 | 0.595723 |
</div>
<br>

A comparison of two different data structures can be seen in <a style="color:inherit;" href="#table2">table 2</a>, here it can be seen that the execution time is significantly faster when using std::vector.

</div>

## 6&emsp;Discussion

#### 6.1&emsp;Ideal gas law

The simulation follows the ideal gas law. As seen in figure 5 the temperature remain almost constant as size increases and pressure drops. With one exception, being that for small sizes, the temperature is lower than it should. This is likley due to that in a smaller space with many particles, there will be a lot more collision, due to our limitations, if a particle collides with another particle, it can not wall bounce in that timestep, resulting in less pressure.

In figure 6, it can be observed that 

#### 6.2&emsp;Superlinear speedup

#### 6.3&emsp;Linear vs Linked List

## 7&emsp;Conclusion

## References

<div id="reference1"></div> 

1. Message Passing Interface :: High Performance computing. (n.d.). https://hpc.nmsu.edu/discovery/mpi/introduction/

</div>