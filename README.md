# GaleShapley Stable Matching Algorithm
In this project, we have a M/M/c queue. The queue entities are considered as 3 different 
surgery requests (A, B and C) joining the queue with rates coming from a poisson distribution.
The servers are 3 different hospital (a, b and c) with service times comming from a normal 
distribution.

Each request ranks the hospitals based on the surgery time (service time of each hospital) and 
each hospital ranks the requests based on the surgery costs (aiming for more profit).

We establish a stable matching based on the hospitals and surgery requests preference relations 
using the [Gale-Shapley Stable Matching Algorithm](https://en.wikipedia.org/wiki/Gale%E2%80%93Shapley_algorithm).
The algorithm is applied on the requests and servers on regular slots (every 8 clocks) and unserved 
requests will get back to the queue waiting for the next slot. The queue is FIFO. We keep the queue 
stable by changing the poisson rates. 

The results include the mean surgery duration for each surgery type, the mean duration of a request
from the time it is issued until the request is completely fullfilled, and the mean cost of each 
surgery type.

Furthermore, we consider a parameter k, which is the number of slots an entity can wait in the 
queue if the best choice for it (the fastest hospital) is not available in that slot. The 
algorithm is run for k = 0, 1,...,5 and the results are plotted. It is shown that as k grows, the 
results improve for both requests and the hospitals.
