The automatic data deleting needs c++-17 so compile with 

g++ -std=c++17 main.cpp -o prog

and run with ./prog

for optimisation try using

g++ -std=c++17 -O3 main.cpp -o prog

or safer

g++ -std=c++17 -O2 main.cpp -o prog

