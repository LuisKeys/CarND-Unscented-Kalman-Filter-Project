#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);
	VectorXd sum(4);
	rmse << 0, 0, 0, 0;

	if (estimations.size() == 0 ||
		estimations.size() != ground_truth.size()) {
		cout << "Incorrect estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		VectorXd estimation = estimations[i];
		VectorXd ground_truth_item = ground_truth[i];

		VectorXd residual = estimation - ground_truth_item;

		VectorXd square = residual.array() * residual.array();
		sum += square;
	}

	//calculate mean
	VectorXd mean = sum / estimations.size();

	//calculate the squared root
	rmse = mean.array().sqrt();

	//return the result
	return rmse;		
}