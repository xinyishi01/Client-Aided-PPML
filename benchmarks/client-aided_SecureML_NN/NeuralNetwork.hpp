#include "defines.hpp"
#include "Layer.hpp"

using namespace std;
using namespace Eigen;

class NeuralNetwork {
public:

	int t1, t2, t3;

	Time_Calculater mytime;

	Mat inputData;
	Mat outputData;
	vector<Layer*> layers;

	NetIO *io;
	int party;
//	ofsteam ofs;

	NeuralNetwork(int _party, NetIO *_io) {
		party = _party;
		io = _io;
		for (int i = 0; i <= LL + 1; ++i) {
			layers.push_back(new Layer(_party));
		}
	}
	~NeuralNetwork() {
		for (vector<Layer*>::iterator it = layers.begin() ; it != layers.end(); ++it)
			delete (*it);

		layers.clear();
	}

	void forward() {
	//	cout << "NN.forward\n";

		(*(layers[0]->getActivation())) = inputData;

	//	cerr << "0\n";
		if (EVAT)
			cerr << "layer: " << 1 << endl;

		if (EVAT)
			t1 = mytime.start();

		layers[1]->forward(inputData);

		if (EVAT)
			cerr << "forward total of L1: " << mytime.end(t1) << endl;

		for (int i = 2; i <= LL; ++i) {
			if (EVAT) {
				t3 = mytime.start();
				cerr << "layer: " << i << endl;
			}
			layers[i]->forward(*(layers[i-1]->getActivation()));
			if (EVAT)
				cerr << "forward total of L" << i << ": " << mytime.end(t3) << endl;
		}

	}

	void backward() {
	//	cout << "NN.backward\n";
		if (EVAT)
			t2 = mytime.start();

		computeGradient();

		if (EVAT)
			cerr << "computeGradient: " << mytime.end(t2) << endl;

		// cout << "computeGradient done" << endl;

		if (EVAT)
			t2 = mytime.start();
		
		updateEquations();

		if (EVAT)
			cerr << "updateEquations: " << mytime.end(t2) << endl;
	}
/*
	myType randomlong(){
		myType rand1 = abs(rand());
	    myType rand2 = abs(rand());
	    rand1 = rand1 << (sizeof(int)*8);   
	    myType randULL = (rand1 | rand2);   
	    return randULL;
	}

	void send_mat(const Mat& x, NetIO *io) {
	    vector<unsigned long int> temp(x.cols()*x.rows());
	    for(int j=0; j<x.rows(); j++) {
	        for(int k=0; k<x.cols(); k++) {
	            temp[j*x.cols()+k] = x(j,k);
	        }
	    }
	    io->send_data(&temp[0], sizeof(unsigned long int) * temp.size());
	    io->flush();
	}

	void receive_mat(Mat& x, NetIO *io) {
	    vector<unsigned long int> temp(x.cols()*x.rows());
	    io->recv_data(&temp[0], sizeof(unsigned long int) * temp.size());
	    for(int j=0; j<x.rows(); j++) {
	        for(int k=0; k<x.cols(); k++) {
	            x(j,k) = temp[j*x.cols()+k];
	        }
	    }
	    //io->flush();
	}*/

	void computeGradient() {
	//	cout << "NN.computeGradient\n";

		if (EVAT)
			t1 = mytime.start();

		int rows = B;
		int columns = LAST_LAYER_SIZE;
		int size = rows * columns;
		int index;

		Mat tmp_ = layers[LL]->y;
		Mat quotient(rows, columns);

	if (QUE) {
		for (int i = 0; i < rows; ++i) {
			Batcher batcher1, batcher2;
			if (party == ALICE) {
				for (int j = 0; j < columns; ++j) {
					batcher1.add<Integer>(64, (unsigned long) tmp_(i, j));
	    			batcher2.add<Integer>(64, (unsigned long) 0);
				}
			}
			else {
				for (int j = 0; j < columns; ++j) {
					batcher1.add<Integer>(64, (unsigned long) 0);
	    			batcher2.add<Integer>(64, (unsigned long) tmp_(i, j));
				}
			}
			batcher1.make_semi_honest(ALICE);
	    	batcher2.make_semi_honest(BOB);
			Integer s = batcher1.next<Integer>();
			s = s + batcher2.next<Integer>();
			for (int j = 1; j < columns; ++j) {
				s = s + batcher1.next<Integer>();
				s = s + batcher2.next<Integer>();
			}

			for (int j = 0; j < columns; ++j) {
				Batcher batcher3, batcher4;
				if (party == ALICE) {
					batcher3.add<Integer>(64, (unsigned long) tmp_(i, j) * (1 << (L)));
	    			batcher4.add<Integer>(64, (unsigned long) 0);
				}
				else {
					batcher3.add<Integer>(64, (unsigned long) 0);
	    			batcher4.add<Integer>(64, (unsigned long) tmp_(i, j) * (1 << (L)));
				}
				batcher3.make_semi_honest(ALICE);
	    		batcher4.make_semi_honest(BOB);
	    		Integer temp1 = batcher3.next<Integer>();
				temp1 = temp1 + batcher4.next<Integer>();
				Integer ans = temp1 / s;
				int pp = 0;
				int p = 1;
				unsigned long x = ans.reveal<long long>(pp);
				Integer r = layers[LL]->rs;
				ans = ans ^ r;
				unsigned long x1 = ans.reveal<long long>(pp);
				unsigned long ra = layers[LL]->r;
				unsigned long ans_ = 0;
				for (int i = 63; i >= 0; --i) {
					unsigned long tmp0 = (party == ALICE ? (((unsigned long)x1 >> i) & 1) : 0) + ((ra >> i) & 1) - 2 * (((unsigned long)x1 >> i) & 1) * ((ra >> i) & 1);
					ans_ = (ans_ * 2) + tmp0;
				}
				quotient(i, j) = (unsigned long)ans_;
			}
		}
	}

		if (EVAT)
			cerr << "division: " << mytime.end(t1) << endl;

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < columns; ++j) {
				(*(layers[LL]->getGradient()))(i, j) = quotient(i, j) - outputData(i, j);
			}

		if (EVAT)
			t1 = mytime.start();

		for (int i = LL; i > 0; --i) {
			if (EVAT)
				cerr << "computeG Layer: " << i << endl;
			layers[i]->computeGradient(*(layers[i-1]->getGradient()));
		}

		if (EVAT)
			cerr << "rest computeGradient: " << mytime.end(t1) << endl;
	}

	void updateEquations() {
	//	cout << "NN.updateEquations\n";

		for (int i = LL; i > 0; --i) {
		//	cout << "----------" << i << "\n";
			layers[i]->updateEquations(*(layers[i-1]->getActivation()));	
		}

	//	layers[0]->updateEquations(inputData);
	}

};
