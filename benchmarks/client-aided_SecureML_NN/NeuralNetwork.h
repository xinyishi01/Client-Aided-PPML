#include "defines.h"
#include "Layer.hpp"

using namespace std;
using namespace Eigen;

class NeuralNetwork {
public:
	Mat inputData;
	Mat outputData;
	vector<Layer*> layers;

	NeuralNetwork(int _party);
	~NeuralNetwork();
	void forward();
	void backward();
	void computeGradient();
	void updateEquations();/*
	void predict(vector<myType> &maxIndex);
	void getAccuracy(const vector<myType> &maxIndex, vector<size_t> &counter);*/
};