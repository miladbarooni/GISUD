#pragma once
#include "ISUD_Base.h"
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include "InitialSolutions.h"
#include "CplexMIP.h"
#include "ISUD.h"


// Construct ISUD_Base from paths
// "mps_file" is the mps_file path
// "rhs_file" is the rhs_file path
// "initial_solution_file" is the initial solution path
// "fixed_cost_file" is the fixed cost file
ISUD_Base generateProblemFromMps(std::string mps_file, std::string rhs_file, std::string initial_solution_file = "", std::string fixed_cost_file = "");
void writeDualsToFile(const std::vector<double>& duals, const std::string& filename);