// Class to Provides various test cases (diffusion, exact solution, and their derivatives)
//
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "TestCase.hpp"
#include "cell.hpp"
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include <Eigen/Dense>

using namespace HArDCore3D;

// ----------------------------------------------------------------------------
//                          Implementation
// ----------------------------------------------------------------------------

// Class
TestCase::TestCase(const std::vector<int> iTC)
	: iTC(iTC),
		_deg_diff(0) {
		validate();
		if (iTC[1]==2){
			_deg_diff = 2;
		}
	}
	
/////////////////// SOLUTION ///////////////////////////////:

// Solution
double TestCase::sol(const double x, const double y, const double z){
	double u = 0;	
	switch(iTC[0]){
			/// iTC[0]=1: \f$u(x,y,z)=sin(\pi x)  sin(\pi y)  sin (\pi z)\f$
		case 1: u = sin(pi*x) * sin(pi*y) * sin(pi*z); 
						break;
			/// iTC[0]=2: \f$u(x,y,z)=cos(\pi x)  cos(\pi y)  cos(\pi z)\f$
		case 2: u = cos(pi*x) * cos(pi*y) * cos(pi*z); 
						break;
			/// iTC[0]=3: \f$u(x,y,z)= x\f$
		case 3: u = x; 
						break;
			/// iTC[0]=4: \f$u(x,y,z)= y\f$
		case 4: u = y; 
						break;
			/// iTC[0]=5: \f$u(x,y,z)= z\f$
		case 5: u = z; 
						break;
		 /// iTC[0]=6: \f$u(x,y,z)= x^2+y^2+z^2\f$
		case 6: u = pow(x,2) + pow(y,2) + pow(z,2);
						break;

		default: break;
	}
	return u;
}

// Gradient of the solution
Eigen::Vector3d TestCase::grad_sol(const double x, const double y, const double z, const Cell* cell){
	Eigen::Vector3d G = Eigen::Vector3d::Zero();
	switch(iTC[0]){
		case 1: G(0) = pi * cos(pi*x) * sin(pi*y) * sin(pi*z);
						G(1) = pi * sin(pi*x) * cos(pi*y) * sin(pi*z); 
						G(2) = pi * sin(pi*x) * sin(pi*y) * cos(pi*z); 
						break;

		case 2: G(0) = -pi * sin(pi*x) * cos(pi*y) * cos(pi*z);
						G(1) = -pi * cos(pi*x) * sin(pi*y) * cos(pi*z);
						G(2) = -pi * cos(pi*x) * cos(pi*y) * sin(pi*z); 
						break;

		case 3: G(0) = 1;
						G(1) = 0;
						G(2) = 0; 
						break;

		case 4: G(0) = 0;
						G(1) = 1;
						G(2) = 0; 
						break;

		case 5: G(0) = 0;
						G(1) = 0;
						G(2) = 1; 
						break;

		case 6: G(0) = 2*x;
						G(1) = 2*y;
						G(2) = 2*z; 
						break;

		default: break;
	}
	return G;
}

// Hessian of the solution
Eigen::Matrix3d TestCase::hess_sol(const double x, const double y, const double z, const Cell* cell){
	Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
	switch(iTC[0]){

		case 1: H.row(0) << - pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z), pi*pi*cos(pi*x)*cos(pi*y)*sin(pi*z), pi*pi*cos(pi*x)*sin(pi*y)*cos(pi*z);
						H.row(1) <<  pi*pi*cos(pi*x)*cos(pi*y)*sin(pi*z), -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z), pi*pi*sin(pi*x)*cos(pi*y)*cos(pi*z);
						H.row(2) << pi*pi*cos(pi*x)*sin(pi*y)*cos(pi*z), pi*pi*sin(pi*x)*cos(pi*y)*cos(pi*z), - pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z);
						break;

		case 2: H.row(0) << - pi*pi*cos(pi*x)*cos(pi*y)*cos(pi*z), pi*pi*sin(pi*x)*sin(pi*y)*cos(pi*z), pi*pi*sin(pi*x)*cos(pi*y)*sin(pi*z);
						H.row(1) <<  pi*pi*sin(pi*x)*sin(pi*y)*cos(pi*z), -pi*pi*cos(pi*x)*cos(pi*y)*cos(pi*z), pi*pi*cos(pi*x)*sin(pi*y)*sin(pi*z);
						H.row(2) << pi*pi*sin(pi*x)*cos(pi*y)*sin(pi*z), pi*pi*cos(pi*x)*sin(pi*y)*sin(pi*z), - pi*pi*cos(pi*x)*cos(pi*y)*cos(pi*z);
						break;

		case 3: break;
		case 4: break;
		case 5: break;

		case 6: H.row(0) << 2, 0, 0;
						H.row(1) << 0, 2, 0;
						H.row(2) << 0, 0, 2;
						break;

		default: break;
	}
	return H;
}

//////////////////////////// DIFFUSION /////////////////////////////

// Diffusion matrix
Eigen::Matrix3d TestCase::diff(const double x, const double y, const double z, const Cell* cell){
	Eigen::Matrix3d K = Eigen::Matrix3d::Identity();
	switch(iTC[1]){
			/// iTC[1]=1: Diff = Id
		case 1: break;		
			/// iTC[1]=2: Diff = \f$\left[\begin{array}{ccc}y^2+z^2+1 & -xy & -xz\\ -xy & x^2+y^2+1 & -yz\\ -xz & -yz  & x^2+y^2+1\end{array}\right]\f$
		case 2: K.row(0) << pow(y,2)+pow(z,2)+1, -x*y , -x*z;			
						K.row(1) << -x*y , pow(x,2)+pow(z,2)+1, -y*z;
						K.row(2) << -x*z, -y*z, pow(x,2)+pow(y,2)+1;
						break;
		default: break;
	}
	return K;
}

// Divergence by row of the diffusion matrix
Eigen::Vector3d TestCase::div_diff(const double x, const double y, const double z, const Cell* cell){
	Eigen::Vector3d divK = Eigen::Vector3d::Zero();
	switch(iTC[1]){
		case 1: break;
		case 2: divK(0) = -2*x;
						divK(1) = -2*y;
						divK(2) = -2*z;
						break;
		default: break;
	}
	return divK;
}

///////////////////////////// SOURCE TERM ///////////////////////////

// Source term
double TestCase::source(const double x, const double y, const double z, const Cell* cell){

	Eigen::Matrix3d AHu = diff(x,y,z,cell) * hess_sol(x,y,z,cell);

	return -AHu.trace() - div_diff(x,y,z,cell).dot(grad_sol(x,y,z,cell));
}



///////////////////////////// VALIDATION ////////////////////////////

void TestCase::validate(){
	
	if (iTC[0]>6 || iTC[1]>3 || (iTC[1]==3 && iTC[0] !=1)){
		std::cout << "Incorrect choice of test cases: iTC= " << iTC[0] << ", " << iTC[1] << "\n";
		exit(EXIT_FAILURE);
	}	

}

