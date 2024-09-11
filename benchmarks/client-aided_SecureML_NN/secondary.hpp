#include "defines.h"

using namespace std;
using namespace Eigen;

NetIO * io;
SHOTIterated* ot;
SHOTIterated* ot2;

Mat train_data(N,D), train_label(N,1);
Mat a[5][5], b[5][5], c[5][5], a_batch[5][5], b_batch[5][5], c_batch[5][5];

void load_train_data(Mat& train_data0, Mat& train_label0, int party) {
    ifstream F( "data"+to_string(party)+".txt" );
    
    cout<<"load training data.......\n";

    int i=0;
    while(i<N) {
        string s;
        if (!getline(F,s)) {
            break;
        }
        istringstream ss(s);
        char c;

        //read data
        for(int j=0; j<D; j++) {
            ss>>train_data0(i,j);
            ss>>c;

        }
        i++;
    }

    i=0;
    while(i<N) {
        string s;
        if (!getline(F,s)) {
            break;
        }
        istringstream ss(s);
        char c;

        //read data
        ss>>train_label0(i);
        ss>>c;

        i++;
    }

    F.close();
    return;
}

void prework(int party) {
    load_train_data(train_data, train_label, party-1);
}

void next_batch(Mat& batch,int start, vector<int>& perm, Mat& A) {
    for (int i = 0; i < B; ++i) {
        batch.row(i) = A.row(perm[start+i]);
    }
    return;
}

void readMiniBatch(NeuralNetwork* net, string phase) {
	next_batch(x_batch, start, perm, train_data);
    next_batch(y_batch, start, perm, train_label);
    net.inputData = x_batch;
    net.outputData = y_batch; // !!!!!!!
    for (int i = 1; i <= LL; ++i)
    	for (int j = 1; j < 5; ++j) {
	    	next_batch(a_batch[i][j], start, perm, a[i][j]);
		    next_batch(b_batch[i][j], start, perm, b[i][j]);
		    next_batch(c_batch[i][j], start, perm, c[i][j]);
    	}
   	for (int i = 1; i <= LL; ++i) {
		net.layers[i]->set_MT(a[i][1], b[i][1], c[i][1], a[i][2], b[i][2], c[i][2], a[i][3], b[i][3], c[i][3], a[i][4], b[i][4], c[i][4]);
	}
}


void train(NeuralNetwork* net, NeuralNetConfig* config) {
	cout << "train\n";,

	for (int i = 0; i < IT; ++i) {
		// cout << "----------------------------------" << endl;  
		// cout << "Iteration " << i << endl;
		
		readMiniBatch(net, "TRAINING");

		// start_m();
		net->forward();
		// end_m("Forward");

		// start_m();
		net->backward();
		// end_m("Backward");

		// cout << "----------------------------------" << endl;  
	}

	// if (STANDALONE)
	// 	net->layers[LL].layerPrint(DEBUG_CONST, "WEIGHTS", DEBUG_PRINT);
	// else
	// 	if (PRIMARY)
	// 		net->layers[LL].funcReconstruct2PC(net->layers[LL].weights, DEBUG_CONST, "weights");
}