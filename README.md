# Evolutionary Structural Optimisation (ESO) algorithm for thermal science

This code is intended to solve the Area-to-point problem in thermal science with an ESO algorithm. This algorithm has never been published in any peer-reviewed journal, since it's very well-known and has been implemented many times before mine. 

It is very easy to use: enter the [filling ratio and the ratio of conductivity of two materials](https://github.com/Raphael-Boichot/Evolutionary-structural-optimisation-algorithm/blob/90a358c182916c89e3a02d2d4ccc7431a01dd4e0/Codes/Run_ESO_method.m#L8) on a heating surface linked to a localized heat sink and it makes the conductive matter (in dark) evolve following a very simple principle :
- find the position of the least quantity (a pixel) of conductive matter that does not impede maximal temperature of the domain when it's removed;
- find the position of the least quantity (a pixel) of conductive matter that decreases at most maximal temperature of the domain when it's added;
- exchange the two positions so that quantity of draining (conductive) material is constant and continue;
- the codes stops when repeated exchanges of positions are detected.
  
The shape obtained presents a very efficient (and optimal) design to cool a distributed heated surfaces like a computer chips, battery stacks, some parts of fuel cells, etc. The theoretical optimal solution to the problem must have equalized temperatures along the adiabatic borders. The code is very simple and super effective, but also super slow to converge. Results are similar to the [Genetic Algorithm case](https://github.com/Raphael-Boichot/A-genetic-algorithm-for-topology-optimization-of-area-to-point-heat-conduction-problem).

**Code free to use, please cite the author according to the license !**

The code is based on an finite difference approximation of the temperature equation solved on arbitrary domains. It uses a sparse direct solver and can converges on moderately powerfull computer within an hour to a day.

## Test case
![](https://github.com/Raphael-Boichot/Evolutionary-structural-optimisation-algorithm/blob/main/Pictures/Test_case.png)

## Exemple of code output

## Exemple of converged shapes
