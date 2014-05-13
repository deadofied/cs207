------
README
------

When operating any of the force-directed models, be sure to zoom in and out and pan the mouse around so that you can see all aspects of convergence. The same applies for the fireworks animation. 

The fireworks animation is just to demonstrate our SDL viewer extensions. Ideally even our force-directed models would look as colorful as the fireworks, but we had some trouble getting colors to work well for these since there were too many nodes/edges/forces for SDL viewer to display colors simultaneously with such a large node size. Sometimes when you zoom in you can see the different colors, e.g. on options (3) and (4).  


--------
DATASETS
--------

Les Miserables: Character co-occurence/coappearence in Victor Hugo's Les Misérables. If character groupings are taken into account by the algorithm then related characters should be placed in closer proximity, while unrelated characters are farther apart, looking something like this: http://bl.ocks.org/mbostock/4062045 :)

US Senators voting: Voting patterns of US senators. Edges connect senators who vote more similarly, based on a correlation threshold of number of votes, i.e. if they've voted (same vote on a bill) together more than a certain number of times, there's an edge.

Karate: Social network of friendships between 34 members of a karate club at a US university in the 1970.



----------
ALGORITHMS
----------


(1) Gauss-Seidel, stabilize in concentric circles (Les Miserables dataset)
## Arrangement without taking chapters/groups into account. We are just looking at the frequency of co-occurence so characters that appear more often (i.e. lots of co-occurences in Le Mis) will have more edges. Gauss-Seidel seeks to bring the arrangement from the initial randomness to more structure / more energetically stable. 

(2) Gauss-Seidel, stabilize in groups (Les Miserables dataset)
## The 7-8 groups are chapters, i.e characters who occurred in the same chapter. 

(3) Fruchterman-Reingold
## 6 nodes were randomly generated. Edges were added between some of them. The progression towards a symmetric arrangement is shown. Relies on spring forces, similar to those in Hooke’s law. There are repulsive forces between all nodes, but also attractive forces between nodes that are adjacent. The label in the SDL viewer is the temperature, which eventually comes down to 0 in a simulated-annealing-esque approach. 

(4) Kamada Kawai with Simulated Annealing Approach
## Instead of a random arrangement, we create 16 nodes that are all connected to their neighbors, so that the ultimate arrangement is a square. Forces between the nodes are computed based on their graph theoretic distances, determined by the lengths of shortest paths between them. The algorithm of Kamada and Kawai uses spring forces proportional to the graph theoretic distances.

(5) Kamada Kawai (Les Miserables dataset)

(6) Gauss-Seidel, stabilize in groups (US Senators voting dataset)
## As you can see, the republicans are on one side and the democrats on another! The two independents are very visible on a separate area of the graph. Ideally this would've looked like this: http://news.yahoo.com/the-splitting-of-the-senate--now-in-convenient-gif-form-213908185.html :)

(7) Kamada Kawai (US Senators voting dataset)
## Shows that Kamada Kawai is not well suited when you are aiming for separation into groups. 

(8) Kamada Kawai (karate dataset)
This nicely arranges itself to show the social network. Demonstrates that Kamada Kawai is ideally suited for social networks.
