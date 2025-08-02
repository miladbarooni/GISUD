// ConsoleApplication1.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include <iostream>
#include "ISUD.h"
#include <string>
#include <ilcplex/ilocplex.h>
#include <cstring>
#include "CplexMIP.h"
#include <Eigen/Dense>
#include "IB_CompatibilityChecker.h"
#include <numeric>
#include "Problems.h"

// Main function
int main(int argc, char* argv[])
{
    std::cout << "Hello" << std::endl;
    std::string path = argv[1];

    bool cplex = strcmp(argv[2], "1")==0;

    ISUD_Base problem = generateProblemFromMps(path + "/columns.txt", path + "/rhs.txt", path + "/initial.txt",
                        path + "/fixed_cost.txt");

    // Solve using GISUD    
    if(! cplex) {

        // Whether to use the column addition strategy in GISUD
        bool colAdd = strcmp(argv[3], "1")==0;
        // Whether to enable the "competition mode". In this mode, 
        // the algorithm follows the addition strategy path but zoom is called
        // each time the column addition strategy
        bool compete = strcmp(argv[4], "1") == 0;

        ISUD isud(&problem, path, colAdd, false, compete);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        std::vector<double> duals;
        isud.solve(duals, path);
        writeDualsToFile(duals, path+"/duals.txt");
        problem.destroyColumns();
    } 
    // Solve using CPLEX
    else {

        ISUD isud(&problem,path, false, false, true);
        std::vector<int> initialSolution = isud.getCurrentSolution();
        CplexMIP cplex(&problem);
        std::vector<int> solution;

        cplex.solve(initialSolution, &solution, false, path);
        problem.destroyColumns();
    }

    return 0;
}
