#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>
#define pi 3.1415926535
#define E  2.71828182845904523536
using namespace std;


void print(vector<double>pos)
{
	for(int i=0;i<pos.size();i++)
	{
		cout << pos[i] << " ";
	}
	cout << endl;

}



double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}

bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}

double fitnessFunction(vector<double> pos)
{
	double result = 0.0;
	for (int i = 0; i < pos.size(); i++)
	{
		result += pow(pos[i], 2);
	}

	return result;
}


class Particle
{
public:

	vector<double> position;

	double VF;
	double F;
	double PF;

	double VCR;
	double CR;
	double PCR;
	bool succ;
	Particle() {}

};

class PLADE
{
public:

	PLADE(int dim, int m, int Tmax, double max, double min, double Fmax, double Fmin,
		double CRmax, double CRmin, double c1, double c2,
		double percent, double wmax, double wmin)
	{

		this->dim = dim;
		this->m = m;
		this->Tmax = Tmax;
		this->max = max;
		this->min = min;
		this->Fmax = Fmax;
		this->Fmin = Fmin;
		this->CRmax = CRmax;
		this->CRmin = CRmin;
		this->c1 = c1;
		this->c2 = c2;
		this->percent = percent;
		this->wmin = wmin;
		this->wmax = wmax;
		particles.resize(m);
		gbest_succ = 0;
		gbest = 0;
		T = 1;
		this->frange = percent * (Fmax - Fmin);
		this->CRrange = percent * (CRmax - CRmin);
	}

	void initialParticles(int i)
	{
		particles[i].position.resize(dim);
		for (int j = 0; j < dim; j++)
		{
			particles[i].position[j] = randDouble(min, max);
			//cout << "particles[i].position:" << particles[i].position[j] << endl;
		}

		particles[i].VF = randDouble(-frange, frange);
		particles[i].VCR = randDouble(-CRrange, CRrange);
		/*cout << "particles[i].VF:" << particles[i].VF << endl;
		cout << "particles[i].VCR:" << particles[i].VCR << endl;*/
		
		particles[i].F = randDouble(Fmin, Fmax);
		particles[i].CR = randDouble(CRmin, CRmax); 
		particles[i].PF = particles[i].F;
		particles[i].PCR = particles[i].CR;
		/*cout << "particles[i].F:" << particles[i].F << endl;
		cout << "particles[i].CR:" << particles[i].CR << endl;*/
	}

	void initialAllParticles()
	{
		for (int i = 0; i < m; i++)
		{
			initialParticles(i);
			if (better(fitnessFunction(particles[i].position), 
				fitnessFunction(particles[gbest].position)))
			{
				gbest = i;
				gbest_succ = i;
			}
		}

	}

	void updateParticle(int i)
	{
		//update f and cr
		
		particles[i].succ = false;
		particles[i].VF = w * particles[i].VF +
			c1 * randDouble(0, 1) * (particles[i].PF - particles[i].F)
			+ c2 * randDouble(0, 1) * (particles[gbest_succ].PF - particles[i].F);

		if (particles[i].VF > frange)
			particles[i].VF = randDouble(-frange, frange);
		else if (particles[i].VF < -frange)
			particles[i].VF = randDouble(-frange, frange);

		particles[i].F = particles[i].F + particles[i].VF;

		if (particles[i].F > Fmax)
			particles[i].F = randDouble(Fmin, Fmax);
		else if (particles[i].F < Fmin)
			particles[i].F = randDouble(Fmin, Fmax);
		/*if (particles[i].F > Fmax)
			particles[i].F = Fmax;
		else if (particles[i].F < Fmin)
			particles[i].F = Fmin;*/

		particles[i].VCR = w * particles[i].VCR +
			c1 * randDouble(0, 1) * (particles[i].PCR - particles[i].CR)
			+ c2 * randDouble(0, 1) * (particles[gbest_succ].PCR - particles[i].CR);

		if (particles[i].VCR > CRrange)
			particles[i].VCR = randDouble(-CRrange, CRrange);
		else if (particles[i].VCR < -CRrange)
			particles[i].VCR = randDouble(-CRrange, CRrange);


		particles[i].CR = particles[i].CR + particles[i].VCR;

		if (particles[i].CR > CRmax)
			particles[i].CR = randDouble(CRmin, CRmax);
		else if (particles[i].CR < CRmin)
			particles[i].CR = randDouble(CRmin, CRmax);
		/*if (particles[i].CR > CRmax)
			particles[i].CR = CRmax;
		else if (particles[i].CR < CRmin)
			particles[i].CR = CRmin;*/

		//update position, in DE processing
		int r1;
		int r2;
		int r3;
		static default_random_engine engine1(time(nullptr));
		static uniform_int_distribution<int> dis1(0, m - 1);

		//mutation
		do
		{
			r1 = dis1(engine1);
			r2 = dis1(engine1);
			r3 = dis1(engine1);

		} while (r1 == r2|| r2 == r3 || r1 == r3 || r1 == i || r2 == i||r3 == i);
		/*cout << "r1:" << r1 << endl;
		cout << "r2:" << r2 << endl;
		cout << "r3:" << r3 << endl;
*/
		vector<double>V = particles[r1].position;
		for (int j = 0; j < dim; j++)
		{
			V[j] = V[j] + particles[i].F*(particles[r2].position[j] - particles[r3].position[j]);
			if (V[j] > max)
				V[j] = max;
			else if (V[j] < min)
				V[j] = min;
		}
		//crossover
		vector<double>U = V;
		static default_random_engine engine2(time(nullptr));
		static uniform_int_distribution<int> dis2(0, dim - 1);
		int rn = dis2(engine2);
		for (int j = 0; j < dim; j++)
		{
			double rd = randDouble(0, 1);
			if (rd > particles[i].CR && j != rn)
			{
				U[j] = particles[i].position[j];
			}
		}

		//seletion
		if (better(fitnessFunction(U), fitnessFunction(particles[i].position)))
		{
			particles[i].position = U;
			particles[i].PCR = particles[i].CR;
			particles[i].PF = particles[i].F;
			particles[i].succ = true;
		}
		if (better(fitnessFunction(particles[i].position), fitnessFunction(particles[gbest].position)))
			gbest = i;

		
	}

	void updateAllParticles()
	{

		for (int i = 0; i < m; i++)
		{
			updateParticle(i);
		}
		
		setBestSucc();
		T++;
	}

	void setBestSucc()
	{
		int index = 0;
		for(int i=0;i<m;i++)
		{
			if(particles[i].succ==true)
			{
				if(better(fitnessFunction(particles[i].position),
					fitnessFunction(particles[index].position)))
				{
					index = i;

				}
			}
		}
		gbest_succ = index;
	}

	void inertiaWeight()
	{
		w = wmax - (wmax - wmin)*T / Tmax;
		cout << "w" << w << endl;
		cout << "particles[i].F:" << particles[gbest_succ].F << endl;
		cout << "particles[i].CR:" << particles[gbest_succ].CR << endl;
		cout << "particles[i].VF:" << particles[gbest_succ].VF << endl;
		cout << "particles[i].VCR:" << particles[gbest_succ].VCR << endl;
		cout << "gbest" << gbest << endl;
		cout << "gbest_succ" << gbest_succ << endl;
		//print(particles[gbest_succ].position);
	}

	int dim;
	int m;
	int Tmax;
	double max;
	double min;

	double w;

	double Fmax;
	double Fmin;

	double CRmax; 
	double CRmin; 
	double wmax;
	double wmin;
	double c1; 
	double c2;
	double percent;
	double frange;
	double CRrange;
	vector<Particle> particles;

	int gbest_succ;
	int gbest;
	int T;

};

void run()
{
	int dim = 30;
	int m = 100;
	int Tmax = 5000;
	double max = 100;
	double min = -100;
	double c1 = 2.0;
	double c2 = 2.0;
	double Fmax = 1;
	double Fmin = 0;
	double CRmax = 1.5;
	double CRmin = -0.5;
	double wmax = 0.9;
	double wmin = 0.4;

	double percent = 0.2;


	PLADE plade = PLADE(dim, m, Tmax, max, min, Fmax, Fmin, CRmax, CRmin, c1, c2,
		percent, wmax, wmin);
	plade.initialAllParticles();
	vector<double>fitness;
	fitness.push_back(fitnessFunction(plade.particles[plade.gbest].position));
	for (int i = 0; i < Tmax; i++)
	{


		plade.inertiaWeight();
		plade.updateAllParticles();
		//plade.caculatePSD();
		//fitness.push_back(pso.getFitness());

		fitness.push_back(fitnessFunction(plade.particles[plade.gbest].position));
		cout << "第" << i << "次迭代结果：";
		cout << ", fitness = " << fitnessFunction(plade.particles[plade.gbest].position) << endl;



	}
}

int main()
{
	run();
	cout << randDouble(-1, 1) << endl;
	system("pause");
}




