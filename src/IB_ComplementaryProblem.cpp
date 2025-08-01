#include "IB_ComplementaryProblem.h"
#include <typeinfo>

// Construct the complementary problem
// "increaseArtificialCost" is a boolean that is true if we want to increase the artificial column cost
// "penalization" is the penalization of the artificialcolumn.
void IB_ComplementaryProblem::constructProblem(bool increaseArtificialCost, double penalization) {
	std::set<int> activeConstraintsIds;
	maxPenalization = penalization;
	int n_tasks = psolutionMethod_->tasks_.size();
	constraints[n_tasks].setBounds(1, 1);
	if (acolId_ != -1) {
		mod.add(constraints[n_tasks]);
	}

	for (int i = 0; i < n_tasks + 1; i++) {
		if (acolId_ == -1) {
			mod.add(constraints[i]);
		}
		if (i < psolutionMethod_->tasks_.size()) {
			constraintsIds[psolutionMethod_->tasks_[i]] = i;
		}
	}

	activeConstraints.add(constraints[n_tasks]);

	double max_p_cost = 0;
	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
		if (psolutionMethod_->columns_[i]->isInCurrentSolution()) {
			if (psolutionMethod_->columns_[i]->getCost() > max_p_cost) {
				max_p_cost = psolutionMethod_->columns_[i]->getCost();
			}
		}
	}


	for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
		IB_Column* column = psolutionMethod_->columns_[i];
		
		if (column->isInCurrentSolution() || (currentPhase_ == -1 || column->getPhase() <= currentPhase_)) {
			int sign = column->isInCurrentSolution() ? -1 : 1;
			double cost_ = column->getCost();
			if (increaseArtificialCost && i == acolId_) {
				cost_ = cost_;
			}
			IloNumColumn col = cost(cost_ * sign);
			for (auto contrib : column->getContribs()) {
				if (contrib.second != 0) {
					col += constraints[constraintsIds[contrib.first]](sign* contrib.second);
					if (!activeConstraintsIds.count(constraintsIds[contrib.first])) {
						activeConstraintsIds.insert(constraintsIds[contrib.first]);
						activeConstraints.add(constraints[constraintsIds[contrib.first]]);
						if (acolId_ != -1) {
							mod.add(constraints[constraintsIds[contrib.first]]);
						}
					}
				}
			}

			 
			if (acolId_ != -1) {

				col += constraints[n_tasks](normalizationConstraint_[i]);
			}
			else {

				col += constraints[n_tasks](normalizationConstraint_[i]);
			}

			colsIndices.push_back(i);
			
			vars.add(IloNumVar(col, 0, acSupport_.count(i) ? 0 : IloInfinity));
			if (i == acolId_) {
				varsacolId = vars.getSize() - 1;
			}

			varsIds[i] = vars.getSize() - 1;
				

			col.end();
		}
	} 

}


// Solve the complementary problem
// Return solution in the vector "solution"
// Past solution can be NULL
// Return dual solution in "duals"
double IB_ComplementaryProblem::solve(std::vector<double>* solution, std::vector<double>* pastSolution, std::vector<double>* duals) {
	if (artificialConstraints.size() > 0) {
		for (int i = 0; i < artificialConstraints.size(); i++) {
			mod.remove(artificialConstraints[i]);
		}
	}

	if (pastSolution != NULL) {
		std::vector<IloExpr> expressions;
		bool first = true;
		IloExpr reference;
		for (int i = 0; i < pastSolution->size(); i++) {
			if ((*pastSolution)[i] > 1e-4) {
				if (first) {
					first = false;
					reference = vars[varsIds[i]];
				}
				else {
					expressions.push_back(reference - vars[varsIds[i]]);
				}
			}
		}
{}{}

		artificialConstraints = std::vector<IloRange>();
		for (int i = 0; i < expressions.size(); i++) {
			artificialConstraints.push_back((expressions[i] == 0));
			mod.add(artificialConstraints[i]);
		}
	}
	
	
	
	if(verbose_level>=2){
		std::cout << "      "<<vars.getSize() << " variables in the CP." << std::endl;
	}
	// Solving the complementary problem
	//
	if (acolId_ != -1) {

		double newObjective = 0;
		double max_value = 0;
		int n = 0;

		IloCplex cplex(mod);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		cplex.setParam(IloCplex::Param::Simplex::Display, 0);
		cplex.setParam(IloCplex::Param::Threads, 8);

		IloNumArray vals(env);
		for (int pen = 0; pen<=5; pen++) {
			cost.setLinearCoef(vars[varsacolId], -(psolutionMethod_->columns_[acolId_]->getCost()) - fabs(psolutionMethod_->columns_[acolId_]->getCost() * 2 * (1./5 * pen)));

			cplex.setOut(env.getNullStream());
			bool success = cplex.solve();
			
			if (success) {
				
				newObjective = 0;
				solution->clear();
				
				cplex.getValues(vars, vals);
				for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
					solution->push_back(0);
				}

				for (int i = 0; i < vals.getSize(); i++) {
					if (vals[i] > 1e-5 && vars[i].getUB() > 0) {
						if (vals[i] > max_value) {
							max_value = vals[i];
						}
						n += 1;
						int sign = psolutionMethod_->columns_[colsIndices[i]]->isInCurrentSolution() ? -1 : 1;
						newObjective += sign * psolutionMethod_->columns_[colsIndices[i]]->getCost() * vals[i];
					}

					solution->at(colsIndices[i]) = vals[i];
				}

				if (vals[varsacolId] > 1e-4 || newObjective >= 0) {
					return newObjective;
				}
				else {

					cplex.setParam(IloCplex::Param::Advance, 1);
					IloNumArray duals(env);
					cplex.getDuals(duals, activeConstraints);
					cplex.setStart(vals, 0, vars, 0, duals, activeConstraints);
				}
			}
			else {
				if(verbose_level>=1){
					std::cout << "      CP Infeasible" << std::endl;
				}
				return 0;
			}
		}

		return newObjective;
	}
	else {
		main_cplex.setParam(IloCplex::Param::Simplex::Display, 0);
		main_cplex.setParam(IloCplex::Param::Threads, 8);
		main_cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-1);
		main_cplex.setOut(env.getNullStream());
		bool success = main_cplex.solve();


		double objective;
		std::set<std::string> coveredTasks;
		if (success) {
			if(verbose_level>=1){
			std::cout << "      CP Success" << std::endl;
			}
			
			IloNumArray vals(env);
			main_cplex.getValues(vars, vals);
			if(duals != NULL) {
				IloNumArray duals_cplex(env);
				main_cplex.getDuals(duals_cplex, constraints);
				for(int i = 0; i < duals_cplex.getSize();i++) {
					duals->push_back(duals_cplex[i]);
				}
			}
			double newObjective = 0;
			int n = 0;
			for (int i = 0; i < psolutionMethod_->columns_.size(); i++) {
				solution->push_back(0);
			}

			double max_value = 0;

			for (int i = 0; i < vals.getSize(); i++) {
				if (vals[i] > 1e-5 && vars[i].getUB() > 0) {
					if (vals[i] > max_value) {
						max_value = vals[i];
					}
					n += 1;
					int sign = psolutionMethod_->columns_[colsIndices[i]]->isInCurrentSolution() ? -1 : 1;
					newObjective += sign * psolutionMethod_->columns_[colsIndices[i]]->getCost() * vals[i];
					for (auto contrib : psolutionMethod_->columns_[colsIndices[i]]->getContribs()) {
						coveredTasks.insert(contrib.first);
					}
				}
				solution->at(colsIndices[i]) = vals[i]; {}
			}

			objective = newObjective;
			int n2 = 0;
			// Suppression des colonnes
			for (int i = 0; i < vars.getSize(); i++) {
				bool delete_ = false;
				for (auto task : coveredTasks) {
					if (psolutionMethod_->columns_[colsIndices[i]]->findContribution(task) > 1e-4) {
						vars[i].setBounds(0, 0);
						delete_ = true;
					}
				}

				if (delete_) {
					n2 += 1;
				}
			}
			if(verbose_level>=2){
				std::cout<<"      " << n2 << " variables put to 0" << std::endl;
			}

		}
		else {
			if(verbose_level>=1){
					std::cout << "      CP Infeasible" << std::endl;
			}
			return 10;
		}

		return objective;
	}
}

// Destroy the complementary problem
void IB_ComplementaryProblem::destroy() {
	env.end();
}
