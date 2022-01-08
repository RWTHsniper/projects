/*
 * simulation.hpp
 *
 *  Created on: Jan 3, 2022
 *      Author: jaeyong
 */

#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <iostream>
#include <Eigen/Dense> // Dense matrices
#include <random> // random sampling
#include<fstream> // write file
#include <memory> // unique_pntr


void buildBrownianMotion(std::vector<Eigen::MatrixXd>& dWIndep, const size_t& numSteps, const size_t& nFactor, 
                                                        const size_t& numPaths, const unsigned int& seed_id = 0, bool normalize = true);
void buildVectofMat(std::vector<Eigen::MatrixXd>& x, const size_t& numSteps, const size_t& nFactor, const size_t& numPaths);
void saveData(std::string fileName, Eigen::MatrixXd  matrix);
std::shared_ptr<Eigen::MatrixXd> computeMSA(std::shared_ptr<Eigen::MatrixXd> r, double dt); // compute money-savings-account


#endif /* SIMULATION_HPP_ */
