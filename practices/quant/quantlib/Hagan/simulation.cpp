#include "simulation.hpp"


void buildBrownianMotion(std::vector<Eigen::MatrixXd>& dWIndep, const size_t& numSteps, const size_t& nFactor, 
                                                        const size_t& numPaths, const unsigned int& seed_id, bool normalize){
    if (numPaths <= 50) normalize = false; // When there are few paths, normalization does not work!
    std::default_random_engine norm_generator{seed_id};
    std::normal_distribution<double> norm_distribution(0.0, 1.0);
    dWIndep.reserve(numSteps);
    Eigen::VectorXd meanSquare(nFactor);
    Eigen::VectorXd variance(nFactor);
    for (size_t i=0; i<numSteps; i++){
        dWIndep.emplace_back(nFactor, numPaths);
        if (normalize) meanSquare.setZero();
        for (size_t j=0; j<nFactor; j++)
            for (size_t k=0; k<numPaths; k++){
                dWIndep[i](j, k) = norm_distribution(norm_generator);
                if (normalize) meanSquare[j] += dWIndep[i](j, k) * dWIndep[i](j, k); // compute the sum of squares
            }
        if (normalize) {
            meanSquare /= numPaths; // compute the mean of squares from sampling
            Eigen::VectorXd meanRow = dWIndep[i].rowwise().mean(); // mean for each state variable
            for (size_t j=0; j<nFactor; j++) variance[j] = meanSquare[j] - meanRow[j] * meanRow[j];
            for (size_t j=0; j<nFactor; j++)
                for (size_t k=0; k<numPaths; k++){
                    dWIndep[i](j, k) = (dWIndep[i](j, k) - meanRow[j]) / std::sqrt(variance[j]);
                }
        }
    }
}

void buildVectofMat(std::vector<Eigen::MatrixXd>& x, const size_t& numSteps, const size_t& nFactor, const size_t& numPaths){
    x.clear();
    x.reserve(numSteps);
    for (size_t i=0; i<numSteps; i++){
        x.emplace_back(nFactor, numPaths);
        x[i].setZero();
    }
}

void saveData(std::string fileName, Eigen::MatrixXd matrix){
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    //https://aleksandarhaber.com/eigen-matrix-library-c-tutorial-saving-and-loading-data-in-from-a-csv-file/
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
 
    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}