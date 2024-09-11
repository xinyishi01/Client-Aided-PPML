//#pragma once
#include "NeuralNetwork.h"
using namespace std;
using namespace Eigen;

NeuralNetwork::NeuralNetwork(int _party) {
	for (int i = 1; i <= LL; ++i) {
		layers.push_back(new Layer(_party));
	}
}


NeuralNetwork::~NeuralNetwork() {
	for (vector<Layer*>::iterator it = layers.begin() ; it != layers.end(); ++it)
		delete (*it);

	layers.clear();
}


void NeuralNetwork::forward() {
	cout << "NN.forward\n";

	layers[1]->forward(inputData);

	for (int i = 2; i <= LL; ++i) {
		layers[i]->forward(*(layers[i-1]->getActivation()));
	}

	// if (PRIMARY)
	// 	funcReconstruct2PC(*(layers[LL]->getActivation()), (*(layers[LL]->getActivation())).size(), "LL: activations");
}

void NeuralNetwork::backward() {
	cout << "NN.backward\n";

	computeGradient();
	// cout << "computeGradient done" << endl;
	updateEquations();
}

void NeuralNetwork::computeGradient() {
	cout << "NN.computeGradient\n";
/*
	int rows = B;
	int columns = LAST_LAYER_SIZE;
	int size = rows*columns;
	int index;

	vector<myType> rowSum(size, 0);
	vector<myType> quotient(size, 0);

	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < columns; ++j)
			rowSum[i*columns] += (*(layers[LL]->getActivation()))(i)(j);

	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < columns; ++j)
			rowSum[i*columns + j] = rowSum[i*columns];

//DIVISION CODE BEGINS HERE
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < columns; ++j) {
			index = i * columns + j;
			if (rowSum[index] != 0)
				quotient[index] = divideMyTypeSA((*(layers[LL]->getActivation()))[index], rowSum[index]);
		}
//DIVISION CODE ENDS HERE

	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < columns; ++j) {
			(*(layers[LL]->getGradient()))(i)(j) = quotient[index] - (outputData(i) == j);
		}

	for (int i = LL; i > 0; --i)
		layers[i]->computeGradient(*(layers[i-1]->getGradient()));*/
}

void NeuralNetwork::updateEquations() {
	cout << "NN.updateEquations\n";

	for (int i = LL; i > 0; --i)
		layers[i]->updateEquations(*(layers[i-1]->getActivation()));	

	layers[0]->updateEquations(inputData);
}