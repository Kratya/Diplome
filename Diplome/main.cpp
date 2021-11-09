#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct points
{
	double x, y, z;
};

const double PI = 3.141593;

class solver
{
public:
	solver() {};
	~solver() {};

	void model_init();
	void model_solve();

private:
	int iter;
	double Func, I, w, A, B, dI, a, b;
//------------------------------------------------------------------
	vector<double> Iq, dIq, I_iq, bq;
	vector<vector<double>> Aqs;
//------------------------------------------------------------------

	int m_iter = 1000;
	double eps = 1e-16, sigma = 0.1, I_i = 1;

	vector<double> V, V_i;
	vector<pair<points, points>> receivers, transmitters;
	double r(points first, points second);
	void init(vector<pair<points, points>>& source, vector<pair<points, points>>& receivers, string fileName);
	double potential(double I, double sigma, pair<points, points> source, pair<points, points> receiver);
	double functional(vector<double> V, vector<double> V_iteration, int size);
	double diff_I(double sigma, pair<points, points> source, pair<points, points> receiver);
	vector<double> calcV(vector<double> Icalc);
	//vector<double> getV();
	vector<double> findDerivativeI(int numRec);
	vector<double> findDerivativeSig();
	void solveI(int indx, vector<double> stval);
	void gSolver(vector<vector<double>> &mat, vector<double> &vec);

	vector<double> calcSmthFunc(vector<double> Tcalc);
	vector<double> findDerivativeT(vector<double> Tcalc, int numRec);
	void solveT(int indx, vector<double> stval);
};

vector<double> solver::calcV(vector<double> Icalc)
{
	vector<double> Vc;
	Vc.resize(receivers.size());

	for (int i = 0; i < receivers.size(); i++)
	{
		points M = receivers[i].first;
		points N = receivers[i].second;
		for (int j = 0; j < transmitters.size(); j++)
		{
			points A = transmitters[j].first;
			points B = transmitters[j].second;
			Vc[i] += (Icalc[j] / (2 * PI * sigma)) * ((1.0 / r(B, M) - 1.0 / r(A, M)) - (1.0 / r(B, N) - 1.0 / r(A, N)));
		}		
	}
	return Vc;
}

vector<double> solver::calcSmthFunc(vector<double> Tcalc)
{
	double A_coeff = 1, B_coeff = 1, C_coeff = 1;
	vector<double> Func;
	Func.resize(receivers.size());

	for (int i = 0; i < receivers.size(); i++)
	{
		points M = receivers[i].first;
		points N = receivers[i].second;
		for (int j = 0; j < transmitters.size(); j++)
		{
			points A = transmitters[j].first;
			points B = transmitters[j].second;
			Func[i] += ((pow(Tcalc[j], 2) * A_coeff + Tcalc[j] * B_coeff + C_coeff)
				/ (2 * PI)) * ((1.0 / r(B, M) - 1.0 / r(A, M)) - (1.0 / r(B, N) - 1.0 / r(A, N)));
		}
	}
	return Func;
}

vector<double> solver::findDerivativeI(int numRec)
{
	vector<double> dVi;
	dVi.resize(receivers.size());

	points A = transmitters[numRec - 1].first;
	points B = transmitters[numRec - 1].second;


	for (int i = 0; i < receivers.size(); i++)
	{
		points M = receivers[i].first;
		points N = receivers[i].second;

		dVi[i] = (1.0 / (2 * PI * sigma)) * ((1.0 / r(B, M) - 1.0 / r(A, M)) - (1.0 / r(B, N) - 1.0 / r(A, N)));
	}
	return dVi;
}

vector<double> solver::findDerivativeSig()
{
	vector<double> dVc;
	dVc.resize(receivers.size());

	for (int i = 0; i < receivers.size(); i++)
	{
		points M = receivers[i].first;
		points N = receivers[i].second;
		for (int j = 0; j < transmitters.size(); j++)
		{
			points A = transmitters[j].first;
			points B = transmitters[j].second;
			dVc[i] += (-Iq[j] / (2 * PI * sigma * sigma)) * ((1.0 / r(B, M) - 1.0 / r(A, M)) - (1.0 / r(B, N) - 1.0 / r(A, N)));
		}
	}
	return dVc;
}

vector<double> solver::findDerivativeT(vector<double> Tcalc, int numRec)
{
	vector<double> dTi;
	double A_coeff = 1, B_coeff = 1;
	dTi.resize(receivers.size());

	points A = transmitters[numRec].first;
	points B = transmitters[numRec].second;


	for (int i = 0; i < receivers.size(); i++)
	{
		points M = receivers[i].first;
		points N = receivers[i].second;

		for (int j = 0; j < transmitters.size(); j++)
		{
			points A = transmitters[j].first;
			points B = transmitters[j].second;

			dTi[i] = ((Tcalc[j] * A_coeff + B_coeff)
				/ (2 * PI)) * ((1.0 / r(B, M) - 1.0 / r(A, M)) - (1.0 / r(B, N) - 1.0 / r(A, N)));
		}

		
	}
	return dTi;
}

void solver::gSolver(vector<vector<double>> &mat, vector<double> &vec)
{
	int n = vec.size();
	double t;

	for (int k = 0; k < n - 1; k++)
	{
		for (int i = k + 1; i < n; i++)
		{
			t = mat[i][k] / mat[k][k];
			vec[i] -= t * vec[k];
			for (int j = k + 1; j < n; j++)
			{
				mat[i][j] -= t * mat[k][j];
			}
		}
	}
	vec[n - 1] /= mat[n - 1][n - 1];

	for (int k = n - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < n; j++)
		{
			sum += mat[k][j] * vec[j];
		}
		vec[k] = (vec[k] - sum) / mat[k][k];
	}
}

void solver::solveT(int indx, vector<double> stval)
{
	vector<double> wq, aVq, last_stval;

	wq.resize(receivers.size());
	V.resize(receivers.size());
	V_i.resize(receivers.size());
	I_iq.resize(transmitters.size());
	last_stval.resize(transmitters.size());
	Aqs.resize(indx);

	double alphaB = 1.0e-14, alpha, FVal = 0, lastFVal, tmp, fdiff;
	vector<vector<double>> matrix;
	matrix.resize(indx);

	for(int i = 0; i < indx; i++)
		matrix[i].resize(indx);

	last_stval = stval;
	V = calcSmthFunc(Iq);
	fill(stval.begin(), stval.end(), 0);

	for (int i = 0; i < indx; i++)
		Aqs[i] = findDerivativeT(Iq, receivers.size());

	for (int i = 0; i < receivers.size(); i++)
		wq[i] = 1.0 / (V[i] * V[i]);

	V_i = calcSmthFunc(I_iq);

	for (int i = 0; i < indx; i++)
		FVal += wq[i] * wq[i] * (V_i[i] - V[i]) * (V_i[i] - V[i]);

	int iter = 0;

	do
	{
		lastFVal = FVal;
		alpha = alphaB;
		V_i = calcSmthFunc(I_iq);

		do
		{
			alpha *= 2;
			for (int i = 0; i < indx; i++)
			{
				for (int j = 0; j < indx; j++)
				{
					matrix[i][j] = 0;
					for (int k = 0; k < indx; k++)
						matrix[i][j] += wq[k] * Aqs[i][k] * Aqs[j][k];
				}
				matrix[i][i] += alpha;
			}
			for (int i = 0; i < indx; i++)
			{
				I_iq = stval;
				fill(stval.begin(), stval.end(), 0);
				for (int k = 0; k < indx; k++)
					stval[i] -= wq[k] * (V_i[i] - V[i]) * Aqs[i][k];
			}
			gSolver(matrix, stval);

			for (int i = 0; i < indx; i++)
				I_iq[i] += stval[i];
			aVq = calcSmthFunc(I_iq);
			for (int i = 0; i < indx; i++)
				I_iq[i] -= stval[i];
			lastFVal = FVal;
			FVal = 0;
			for (int i = 0; i < indx; i++)
				FVal += wq[i] * wq[i] * (aVq[i] - V[i]) * (aVq[i] - V[i]);
		} while (FVal < lastFVal);
		FVal = lastFVal;

		for (int i = 0; i < indx; i++)
			I_iq[i] += last_stval[i];
		tmp = 0;
		for (int i = 0; i < indx; i++)
			tmp += last_stval[i] * last_stval[i];
		fdiff = abs(FVal - lastFVal);
		iter++;
		// Write smth
	} while (tmp > eps * eps || fdiff > eps);

	for (int i = 0; i < indx; i++)
		cout << "F[" << i << "] = " << I_iq[i] << endl;
}

void solver::solveI(int indx, vector<double> stval)
{
	vector<double> wq, aVq, last_stval;

	wq.resize(receivers.size());
	V.resize(receivers.size());
	V_i.resize(receivers.size());
	I_iq.resize(transmitters.size());
	last_stval.resize(transmitters.size());
	Aqs.resize(indx);

	double alphaB = 1.0e-14, alpha, FVal = 0, lastFVal, tmp, fdiff;
	vector<vector<double>> matrix;
	matrix.resize(indx);

	for (int i = 0; i < indx; i++)
		matrix[i].resize(indx);

	last_stval = stval;
	V = calcV(Iq);
	fill(stval.begin(), stval.end(), 0);

	for (int i = 0; i < indx; i++)
		Aqs[i] = findDerivativeI(receivers.size());

	for (int i = 0; i < receivers.size(); i++)
		wq[i] = 1.0 / (V[i] * V[i]);

	V_i = calcV(I_iq);

	for (int i = 0; i < indx; i++)
		FVal += wq[i] * wq[i] * (V_i[i] - V[i]) * (V_i[i] - V[i]);

	int iter = 0;

	do
	{
		lastFVal = FVal;
		alpha = alphaB;
		V_i = calcV(I_iq);

		do
		{
			alpha *= 2;
			for (int i = 0; i < indx; i++)
			{
				for (int j = 0; j < indx; j++)
				{
					matrix[i][j] = 0;
					for (int k = 0; k < indx; k++)
						matrix[i][j] += wq[k] * Aqs[i][k] * Aqs[j][k];
				}
				matrix[i][i] += alpha;
			}
			for (int i = 0; i < indx; i++)
			{
				I_iq = stval;
				fill(stval.begin(), stval.end(), 0);
				for (int k = 0; k < indx; k++)
					stval[i] -= wq[k] * (V_i[i] - V[i]) * Aqs[i][k];
			}
			gSolver(matrix, stval);

			for (int i = 0; i < indx; i++)
				I_iq[i] += stval[i];
			aVq = calcV(I_iq);
			for (int i = 0; i < indx; i++)
				I_iq[i] -= stval[i];
			lastFVal = FVal;
			FVal = 0;
			for (int i = 0; i < indx; i++)
				FVal += wq[i] * wq[i] * (aVq[i] - V[i]) * (aVq[i] - V[i]);
		} while (FVal < lastFVal);
		FVal = lastFVal;

		for (int i = 0; i < indx; i++)
			I_iq[i] += last_stval[i];
		tmp = 0;
		for (int i = 0; i < indx; i++)
			tmp += last_stval[i] * last_stval[i];
		fdiff = abs(FVal - lastFVal);
		iter++;
		// Write smth
	} while (tmp > eps * eps || fdiff > eps);

	for (int i = 0; i < indx; i++)
		cout << "I[" << i << "] = " << I_iq[i] << endl;
}

double solver::r(points first, points second)
{
	return sqrt
	(
		(first.x - second.x) * (first.x - second.x) +
		(first.y - second.y) * (first.y - second.y) +
		(first.z - second.z) * (first.z - second.z)
	);
}

void solver::init(vector <pair<points, points>>& source, vector<pair<points, points>>& receivers, string fileName)
{
	int count;
	points A, B, M, N;
	double I;

	ifstream file(fileName);

	file >> count;

	for (int i = 0; i < count; i++)
	{
		file >> M.x >> M.y >> M.z >> N.x >> N.y >> N.z;
		receivers.push_back(pair<points, points>(M, N));
	}

	file >> count;
	for (int i = 0; i < count; i++)
	{
		file >> A.x >> A.y >> A.z >> B.x >> B.y >> B.z;
		source.push_back(pair<points, points>(A, B));
		file >> I;
		Iq.push_back(I);
	}
	file.close();
}

void solver::model_init()
{
	init(transmitters, receivers, "param.txt");
}

void solver::model_solve()
{
	solveI(Iq.size(), Iq);
}

void main()
{
	solver sl;
	sl.model_init();
	sl.model_solve();

	system("pause");
}
