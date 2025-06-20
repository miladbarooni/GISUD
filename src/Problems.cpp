#include "Problems.h"

// "rhs_file" is the rhs_file path
// "initial_solution_file" is the initial solution path
// "fixed_cost_file" is the fixed cost file
ISUD_Base generateProblemFromMps(std::string mps_file, std::string rhs_file, std::string initial_solution_file, std::string fixed_cost_file) {
    std::ifstream file, rfile, isfile;
    if (rhs_file != "") {
        rfile.open(rhs_file);
    }
    std::cout << mps_file << std::endl;
    file.open(mps_file);
    std::string nRowsS, nColsS;
    file >> nRowsS;
    file >> nColsS;
    std::cout << nRowsS << " " << nColsS << std::endl;
    int nRows = std::stoi(nRowsS);
    std::vector<std::string> tasks;
    std::vector<int> rhs;
    for (int i = 0; i < nRows; i++) {
        tasks.push_back("task_" + std::to_string(i));
    }

    for (int i = 0; i < nRows; i++) {
        std::string rhs_;

        if (rhs_file != "") {
            rfile >> rhs_;
        }
        else {
            rhs_ = "1";
        }

        rhs.push_back(stoi(rhs_));
    }

    std::cout << rhs.size() << " " << tasks.size() << std::endl;
    std::vector<IB_Column*> columns;
    Eigen::MatrixXi matrix(nRows, std::stoi(nColsS));
    matrix = Eigen::MatrixXi::Zero(nRows, std::stoi(nColsS));
    Eigen::VectorXi rhs_vector = Eigen::Map<Eigen::VectorXi>(rhs.data(), rhs.size());
    double max_cost = 0;
    while (file.good() && columns.size() < std::stoi(nColsS)) {

        std::string costS;
        std::string nNonZerosS;

        file >> costS;
        file >> nNonZerosS;
        std::string colName = "column_" + std::to_string(columns.size());

        double cost = std::stod(costS);
        std::vector<std::pair<std::string, int>> contribs;
        std::set<std::string> uniqueContribs;
        for (int i = 0; i < std::stoi(nNonZerosS); i++) {
            std::string rowContrib;
            file >> rowContrib;

            if (!uniqueContribs.count(tasks[((int)stoi(rowContrib))])) {
                contribs.push_back(std::pair<std::string, int>(tasks[((int)stoi(rowContrib))], 1));
                uniqueContribs.insert(tasks[((int)stoi(rowContrib))]);
                matrix(((int)stoi(rowContrib)), columns.size()) = 1;
            }
        }


        columns.push_back(new IB_Column(colName, contribs));
        if (cost > max_cost) {
            max_cost = cost;
        }
        columns[columns.size() - 1]->setCost(cost);
        columns[columns.size() - 1]->setOutCurrentSolution();
        if (initial_solution_file == "") {
            if (contribs.size() > 0 && ((rhs_vector - matrix(Eigen::seq(0, Eigen::last), columns.size() - 1)).array() >= 0).all()) {
                rhs_vector -= matrix(Eigen::seq(0, Eigen::last), columns.size() - 1);
                columns[columns.size() - 1]->setInCurrentSolution();
            }
            else {
                columns[columns.size() - 1]->setOutCurrentSolution();
            }
        }

    }

    file.close();

    int n_artificial_columns = 0;

    std::vector<int> remaining_rhs;
    for (int i = 0; i < rhs_vector.size(); i++) {
        remaining_rhs.push_back(rhs_vector(i));
    }

    std::cout << "Computation of artificial columns." << std::endl;

    if (initial_solution_file == "") {
        if (true) {
            for (auto column : addArtificialColumns(tasks, remaining_rhs, max_cost)) {
                columns.push_back(column);
                n_artificial_columns++;
            }
        }

    }
    else {
        isfile.open(initial_solution_file);
        int n_pos_cols = 0;
        std::set<int> column_pos_nos;
        std::map<int, int> duplicates;
        while (isfile.good()) {
            std::string column_no;
            isfile >> column_no;
            columns[stoi(column_no)]->setInCurrentSolution();
            if (column_pos_nos.count(stoi(column_no))) {
                columns.push_back(new IB_Column(*columns[stoi(column_no)]));
                columns[columns.size() - 1]->setName(columns[columns.size() - 1]->getName() + "_DUPLICATE_" + std::to_string(duplicates[stoi(column_no)]));
                duplicates[stoi(column_no)] += 1;
                columns[columns.size() - 1]->setInCurrentSolution();
            }
            else {
                column_pos_nos.insert(stoi(column_no));
                duplicates[stoi(column_no)] = 1;
            }

            rhs_vector -= matrix(Eigen::seq(0, Eigen::last), stoi(column_no));
            n_pos_cols += 1;
        }

        std::cout << n_pos_cols << " positive columns." << std::endl;


        std::vector<std::string> finalTasks;
        std::vector<int> finalRhs;
        for (int i = 0; i < rhs_vector.size(); i++) {
            if (rhs_vector(i) == 0) {
                finalTasks.push_back(tasks[i]);
                finalRhs.push_back(rhs[i]);
            }
            else {
                std::cout << "Deletion of task : " << tasks[i] << std::endl;
                for (auto column : columns) {
                    column->deleteContrib(tasks[i]);
                }
            }
        }

        tasks = finalTasks;
        rhs = finalRhs;
    }


    //std::cout << rhs_vector.transpose() << std::endl;

    std::cout << n_artificial_columns << " artificial columns." << std::endl;

    std::vector<IB_Column*> final_columns;
    for (auto column : columns) {
        if (column->getContribs().size() > 0 || column->isInCurrentSolution()) {
            column->reorderContribs(tasks);
            final_columns.push_back(column);
        }
    }

    double fixed_cost = 0;
    if (fixed_cost_file != "") {
        std::ifstream fixed_cost_f;
        fixed_cost_f.open(fixed_cost_file);
        std::string fixed_cost_s;
        fixed_cost_f >> fixed_cost_s;
        fixed_cost = stod(fixed_cost_s);
        fixed_cost_f.close();
        std::cout << "Fixed cost : " << fixed_cost << std::endl;
    }

    return ISUD_Base(tasks, rhs, final_columns, fixed_cost);
}