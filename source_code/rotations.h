#ifndef ROTATIONS_H
#define ROTATIONS_H

#include <cmath>
#include <vector>

class ZXZEulerRM
{
public:
	double R[3][3];
	double RT[3][3];

	ZXZEulerRM(const double& alpha, const double& beta, const double& gamma)
	{
		double c1 = cos(alpha);
		double s1 = sin(alpha);
		double c2 = cos(beta);
		double s2 = sin(beta);
		double c3 = cos(gamma);
		double s3 = sin(gamma);
		// Rotation matrix (active rotation)
		R[0][0] = c1*c3 - c2*s1*s3;
		R[0][1] = -c1*s3 - c2*c3*s1;
		R[0][2] = s1*s2;
		R[1][0] = c3*s1 + c1*c2*s3;
		R[1][1] = c1*c2*c3 - s1*s3;
		R[1][2] = -c1*s2;
		R[2][0] = s2*s3;
		R[2][1] = c3*s2;
		R[2][2] = c2;
		// Trnasposed rotation matrix (passive rotation)
		RT[0][0] = R[0][0];
		RT[0][1] = R[1][0];
		RT[0][2] = R[2][0];
		RT[1][0] = R[0][1];
		RT[1][1] = R[1][1];
		RT[1][2] = R[2][1];
		RT[2][0] = R[0][2];
		RT[2][1] = R[1][2];
		RT[2][2] = R[2][2];
	}

	std::vector<double> RxV(std::vector<double> const& v, bool const& t) const 
	{
		std::vector<double> p; p.reserve(3);
		for (size_t i = 0; i < 3; i++) {
			if (t) {
				p.push_back(RT[i][0] * v[0] + RT[i][1] * v[1] + RT[i][2] * v[2]);
			}
			else {
				p.push_back(R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2]);
			}
		}
		return p;
	}

	void RxM(double const M[3][3], double P[3][3], bool const& t) const
	{
		for (size_t i = 0; i < 3; ++i) {
			for (size_t j = 0; j < 3; ++j) {
				if (t) {
					P[i][j] = RT[i][0] * M[0][j] + RT[i][1] * M[1][j] + RT[i][2] * M[2][j];
				}
				else {
					P[i][j] = R[i][0] * M[0][j] + R[i][1] * M[1][j] + R[i][2] * M[2][j];
				}
			}
		}
	}
};

#endif